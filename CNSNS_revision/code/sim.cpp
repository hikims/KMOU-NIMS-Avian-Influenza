// =====================================================================
//  sim.cpp  --  revised simulator for the confined SI reaction-diffusion
//               model (CNSNS revision, October 2026)
//
//  Two transmission models are supported:
//
//    fixed  :  F(x,t) = k1 * I(x,t) * chi(x; r0)      (localized-source model,
//              chi centred at the fixed seed location x0 = domain centre)
//
//    local  :  F(x,t) = k1 * I(x,t)                   (local mass-action model,
//              the delta-kernel limit of the kernel model)
//
//    kernel :  F(x,t) = k1 * (K * I)(x,t)             (nonlocal kernel model,
//              K(r) = chi(r)/Z normalised to unit mass, convolution with
//              reflective (Neumann-consistent) padding evaluated by FFT)
//
//  with  S_t = -S F ,  I_t = D Laplace(I) + S F ,  zero-flux boundary.
//
//  Time stepping: Peaceman-Rachford ADI (Crank-Nicolson accuracy) for the
//  diffusion of I, explicit reaction, S updated with HALF steps in each
//  ADI sweep so that S + I is conserved (the submitted code applied a full
//  dt to S in each half sweep, i.e. S decayed at twice the stated rate).
//
//  Build:  g++ -std=c++17 -O3 -I/opt/homebrew/include -L/opt/homebrew/lib
//              sim.cpp -ltbb -o sim
//  Usage:  ./sim --model fixed --L 15 --D 0.01 --out OUTDIR [options]
//  Options (defaults):
//     --dx 0.04  --dt 0.001  --T 400  --r0 0.5  --eps 1.0  --k1 1
//     --init_i 0.1  --leng 0.25   --thetas 0.005,0.01,0.02,0.05,0.1,0.2,0.3
//     --snap t1,t2,...   --out_every 0.2   --conv_every 10   --threads 0
// =====================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>

typedef std::complex<double> cplx;

// ---------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------
static void mkdir_p(const std::string& path) {
    std::string cur;
    for (size_t i = 0; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            if (!cur.empty()) mkdir(cur.c_str(), 0755);
        }
        if (i < path.size()) cur += path[i];
    }
}

static std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) if (!tok.empty()) v.push_back(atof(tok.c_str()));
    return v;
}

// ---------------------------------------------------------------------
// radix-2 complex FFT (in place, length n = power of two)
// ---------------------------------------------------------------------
struct FFT {
    int n, logn;
    std::vector<int> rev;
    std::vector<cplx> w;   // twiddles for forward transform
    explicit FFT(int n_) : n(n_) {
        logn = 0; while ((1 << logn) < n) ++logn;
        rev.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            int r = 0;
            for (int b = 0; b < logn; ++b) if (i & (1 << b)) r |= 1 << (logn - 1 - b);
            rev[i] = r;
        }
        w.assign(n / 2, cplx(0, 0));
        for (int k = 0; k < n / 2; ++k) w[k] = std::polar(1.0, -2.0 * M_PI * k / n);
    }
    void run(cplx* a, int stride, cplx* buf, bool inverse) const {
        for (int i = 0; i < n; ++i) buf[i] = a[rev[i] * stride];
        for (int len = 2; len <= n; len <<= 1) {
            int half = len >> 1, step = n / len;
            for (int i = 0; i < n; i += len) {
                for (int j = 0; j < half; ++j) {
                    cplx tw = w[j * step];
                    if (inverse) tw = std::conj(tw);
                    cplx u = buf[i + j], v = buf[i + j + half] * tw;
                    buf[i + j] = u + v;
                    buf[i + j + half] = u - v;
                }
            }
        }
        if (inverse) { double s = 1.0 / n; for (int i = 0; i < n; ++i) buf[i] *= s; }
        for (int i = 0; i < n; ++i) a[i * stride] = buf[i];
    }
};

static void fft2d(std::vector<cplx>& A, int P, const FFT& fft, bool inverse) {
    // rows
    tbb::parallel_for(tbb::blocked_range<int>(0, P), [&](const tbb::blocked_range<int>& r) {
        std::vector<cplx> buf(P);
        for (int i = r.begin(); i < r.end(); ++i) fft.run(&A[(size_t)i * P], 1, buf.data(), inverse);
    });
    // columns
    tbb::parallel_for(tbb::blocked_range<int>(0, P), [&](const tbb::blocked_range<int>& r) {
        std::vector<cplx> buf(P);
        for (int j = r.begin(); j < r.end(); ++j) fft.run(&A[j], P, buf.data(), inverse);
    });
}

// ---------------------------------------------------------------------
// tridiagonal (Thomas) with precomputed coefficients
//   a: sub, b: diag, c: super (constant along the line)
// ---------------------------------------------------------------------
struct Tridiag {
    int n;
    std::vector<double> cp, m;    // modified super-diagonal, pivots
    std::vector<double> a;        // sub-diagonal
    Tridiag() : n(0) {}
    void setup(int n_, double r) {
        // (1+2r) x_i - r x_{i-1} - r x_{i+1} = d_i ; Neumann: (1+r) at ends
        n = n_;
        std::vector<double> b(n, 1.0 + 2.0 * r), c(n, -r);
        a.assign(n, -r);
        b[0] = 1.0 + r; b[n - 1] = 1.0 + r;
        cp.assign(n, 0.0); m.assign(n, 0.0);
        m[0] = b[0]; cp[0] = c[0] / m[0];
        for (int i = 1; i < n; ++i) {
            m[i] = b[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / m[i];
        }
    }
    // solve in place: d (stride) -> x
    void solve(double* d, int stride, double* tmp) const {
        tmp[0] = d[0] / m[0];
        for (int i = 1; i < n; ++i) tmp[i] = (d[i * stride] - a[i] * tmp[i - 1]) / m[i];
        d[(size_t)(n - 1) * stride] = tmp[n - 1];
        for (int i = n - 2; i >= 0; --i) d[(size_t)i * stride] = tmp[i] - cp[i] * d[(size_t)(i + 1) * stride];
    }
};

// ---------------------------------------------------------------------
int main(int argc, char** argv) {
    std::string model = "fixed", out = "out";
    double L = 15.0, D = 0.01, dx = 0.04, dt = 0.001, T = 400.0;
    double r0 = 0.5, eps = 1.0, k1 = 1.0, init_i = 0.1, leng = 0.25;
    double out_every = 0.2;
    int conv_every = 10, threads = 0;
    bool consistent_bc = false;
    std::vector<double> thetas = {0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3};
    std::vector<double> snaps;

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&]() -> std::string { if (i + 1 < argc) return argv[++i]; fprintf(stderr, "missing value for %s\n", k.c_str()); exit(1); };
        if (k == "--model") model = next();
        else if (k == "--out") out = next();
        else if (k == "--L") L = atof(next().c_str());
        else if (k == "--D") D = atof(next().c_str());
        else if (k == "--dx") dx = atof(next().c_str());
        else if (k == "--dt") dt = atof(next().c_str());
        else if (k == "--T") T = atof(next().c_str());
        else if (k == "--r0") r0 = atof(next().c_str());
        else if (k == "--eps") eps = atof(next().c_str());
        else if (k == "--k1") k1 = atof(next().c_str());
        else if (k == "--init_i") init_i = atof(next().c_str());
        else if (k == "--leng") leng = atof(next().c_str());
        else if (k == "--out_every") out_every = atof(next().c_str());
        else if (k == "--conv_every") conv_every = atoi(next().c_str());
        else if (k == "--threads") threads = atoi(next().c_str());
        else if (k == "--consistent_bc") consistent_bc = true;
        else if (k == "--thetas") thetas = parse_list(next());
        else if (k == "--snap") snaps = parse_list(next());
        else { fprintf(stderr, "unknown option %s\n", k.c_str()); return 1; }
    }
    if (model != "fixed" && model != "kernel" && model != "local") { fprintf(stderr, "model must be fixed|kernel|local\n"); return 1; }

    std::unique_ptr<tbb::global_control> gc;
    if (threads > 0) gc.reset(new tbb::global_control(tbb::global_control::max_allowed_parallelism, threads));

    const int Nx = (int)llround(L / dx), Ny = Nx;
    const double dy = dx;
    const double area_domain = L * L;
    const size_t NN = (size_t)Nx * Ny;
    const long nsteps = (long)llround(T / dt);
    const int out_stride = std::max(1, (int)llround(out_every / dt));
    const double cx = 0.5 * L, cy = 0.5 * L;

    mkdir_p(out);

    std::vector<double> S(NN), I(NN), Sn(NN), In(NN), F(NN), Fn(NN), G(NN, 0.0);

    // initial condition: disc of radius leng with I = init_i, S = 1 - I
    for (int i = 0; i < Nx; ++i) for (int j = 0; j < Ny; ++j) {
        double x = i * dx, y = j * dy;
        double d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
        double inf = (d2 <= leng * leng) ? init_i : 0.0;
        I[i * Ny + j] = inf; S[i * Ny + j] = 1.0 - inf;
    }

    // fixed cut-off profile chi(x) (used by the fixed model)
    std::vector<double> chi(NN, 1.0);
    for (int i = 0; i < Nx; ++i) for (int j = 0; j < Ny; ++j) {
        double x = i * dx, y = j * dy;
        double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
        chi[i * Ny + j] = 0.5 * (1.0 + std::tanh((r0 - r) / eps));
    }

    // kernel set-up (kernel model)
    int Rk = 0, P = 1;
    std::vector<cplx> Khat, Abuf;
    std::unique_ptr<FFT> fft;
    if (model == "kernel") {
        double Rphys = r0 + 4.0 * eps;             // truncation radius
        Rk = (int)std::ceil(Rphys / dx);
        int need = Nx + 2 * Rk;
        while (P < need) P <<= 1;
        fft.reset(new FFT(P));
        // kernel on the grid, normalised to unit mass
        std::vector<double> ker((size_t)(2 * Rk + 1) * (2 * Rk + 1));
        double Z = 0.0;
        for (int a = -Rk; a <= Rk; ++a) for (int b = -Rk; b <= Rk; ++b) {
            double r = std::sqrt((double)(a * a + b * b)) * dx;
            double v = (r <= Rphys) ? 0.5 * (1.0 + std::tanh((r0 - r) / eps)) : 0.0;
            ker[(size_t)(a + Rk) * (2 * Rk + 1) + (b + Rk)] = v;
            Z += v;
        }
        Z *= dx * dy;
        Khat.assign((size_t)P * P, cplx(0, 0));
        for (int a = -Rk; a <= Rk; ++a) for (int b = -Rk; b <= Rk; ++b) {
            int ia = (a + P) % P, ib = (b + P) % P;
            Khat[(size_t)ia * P + ib] = ker[(size_t)(a + Rk) * (2 * Rk + 1) + (b + Rk)] / Z * dx * dy;
        }
        fft2d(Khat, P, *fft, false);
        Abuf.assign((size_t)P * P, cplx(0, 0));
        printf("kernel: Rk=%d cells, P=%d, Z=%.6f\n", Rk, P, Z);
    }

    // reflective index helper
    auto refl = [](int idx, int n) { if (idx < 0) idx = -idx - 1; if (idx >= n) idx = 2 * n - 1 - idx; return idx; };

    auto convolve = [&](const std::vector<double>& Iin, std::vector<double>& Gout) {
        // build reflect-padded array into Abuf
        tbb::parallel_for(tbb::blocked_range<int>(0, P), [&](const tbb::blocked_range<int>& r) {
            for (int a = r.begin(); a < r.end(); ++a) {
                int i = a - Rk;
                bool inside_i = (i >= -Rk && i < Nx + Rk);
                int ii = inside_i ? refl(i, Nx) : -1;
                for (int b = 0; b < P; ++b) {
                    int j = b - Rk;
                    bool inside = inside_i && (j >= -Rk && j < Ny + Rk);
                    Abuf[(size_t)a * P + b] = inside ? cplx(Iin[(size_t)ii * Ny + refl(j, Ny)], 0.0) : cplx(0, 0);
                }
            }
        });
        fft2d(Abuf, P, *fft, false);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, (size_t)P * P), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t q = r.begin(); q < r.end(); ++q) Abuf[q] *= Khat[q];
        });
        fft2d(Abuf, P, *fft, true);
        tbb::parallel_for(tbb::blocked_range<int>(0, Nx), [&](const tbb::blocked_range<int>& r) {
            for (int i = r.begin(); i < r.end(); ++i)
                for (int j = 0; j < Ny; ++j)
                    Gout[(size_t)i * Ny + j] = std::max(0.0, Abuf[(size_t)(i + Rk) * P + (j + Rk)].real());
        });
    };

    auto compute_F = [&](const std::vector<double>& Iin, std::vector<double>& Fout) {
        if (model == "fixed") {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, NN), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t q = r.begin(); q < r.end(); ++q) Fout[q] = k1 * Iin[q] * chi[q];
            });
        } else if (model == "local") {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, NN), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t q = r.begin(); q < r.end(); ++q) Fout[q] = k1 * Iin[q];
            });
        } else {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, NN), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t q = r.begin(); q < r.end(); ++q) Fout[q] = k1 * G[q];
            });
        }
    };

    const double rx = 0.5 * D * dt / (dx * dx), ry = 0.5 * D * dt / (dy * dy);
    Tridiag tx, ty; tx.setup(Nx, rx); ty.setup(Ny, ry);

    // ---- output files ----
    std::string ts_file = out + "/timeseries.txt";
    FILE* fts = fopen(ts_file.c_str(), "w");
    fprintf(fts, "# model=%s L=%g D=%g dx=%g dt=%g T=%g Nx=%d r0=%g eps=%g k1=%g init_i=%g leng=%g\n",
            model.c_str(), L, D, dx, dt, T, Nx, r0, eps, k1, init_i, leng);
    fprintf(fts, "# thetas:"); for (double th : thetas) fprintf(fts, " %g", th); fprintf(fts, "\n");
    fprintf(fts, "# t massI massS maxI centerI PR var");
    for (double th : thetas) fprintf(fts, " area(%g)", th);
    fprintf(fts, "\n");
    fclose(fts);

    size_t nth = thetas.size();
    std::vector<double> peak_area(nth, -1.0), peak_time(nth, 0.0), final_area(nth, 0.0);
    double massI = 0, massS = 0, maxI = 0, PR = 0, var = 0;
    size_t snap_idx = 0;
    std::sort(snaps.begin(), snaps.end());

    auto write_snapshot = [&](double t) {
        char name[512];
        snprintf(name, sizeof(name), "%s/snap_t%08.2f.bin", out.c_str(), t);
        FILE* fs = fopen(name, "wb");
        fwrite(I.data(), sizeof(double), NN, fs);
        fclose(fs);
        snprintf(name, sizeof(name), "%s/snapS_t%08.2f.bin", out.c_str(), t);
        fs = fopen(name, "wb");
        fwrite(S.data(), sizeof(double), NN, fs);
        fclose(fs);
    };

    auto observe = [&](long k) {
        double t = k * dt;
        double mI = 0, mS = 0, mx = 0, s2 = 0;
        std::vector<long> cnt(nth, 0);
        for (size_t q = 0; q < NN; ++q) {
            double v = I[q];
            mI += v; mS += S[q]; s2 += v * v; if (v > mx) mx = v;
            for (size_t a = 0; a < nth; ++a) if (v >= thetas[a]) ++cnt[a];
        }
        mI *= dx * dy; mS *= dx * dy; s2 *= dx * dy;
        massI = mI; massS = mS; maxI = mx;
        PR = (s2 > 0) ? (mI * mI) / (area_domain * s2) : 0.0;
        var = s2 / area_domain - (mI / area_domain) * (mI / area_domain);
        FILE* f = fopen(ts_file.c_str(), "a");
        fprintf(f, "%.6f %.10f %.10f %.10f %.10f %.10f %.12e", t, mI, mS, mx, I[(size_t)(Nx / 2) * Ny + Ny / 2], PR, var);
        for (size_t a = 0; a < nth; ++a) {
            double ar = cnt[a] * dx * dy;
            fprintf(f, " %.6f", ar);
            if (ar > peak_area[a]) { peak_area[a] = ar; peak_time[a] = t; }
            final_area[a] = ar;
        }
        fprintf(f, "\n");
        fclose(f);
        while (snap_idx < snaps.size() && snaps[snap_idx] <= t + 0.5 * dt) { write_snapshot(t); ++snap_idx; }
    };

    printf("model=%s L=%g D=%g Nx=%d dt=%g steps=%ld\n", model.c_str(), L, D, Nx, dt, nsteps);
    if (model == "kernel") convolve(I, G);
    observe(0);

    for (long k = 1; k <= nsteps; ++k) {
        if (model == "kernel" && ((k - 1) % conv_every == 0) && k > 1) convolve(I, G);
        compute_F(I, F);

        // ---- half step 1: explicit in y, implicit in x ----
        tbb::parallel_for(tbb::blocked_range<int>(0, Ny), [&](const tbb::blocked_range<int>& r) {
            std::vector<double> tmp(Nx);
            for (int j = r.begin(); j < r.end(); ++j) {
                int jp = (j == Ny - 1) ? (consistent_bc ? j : j - 1) : j + 1;
                int jm = (j == 0) ? (consistent_bc ? j : j + 1) : j - 1;
                for (int i = 0; i < Nx; ++i) {
                    size_t q = (size_t)i * Ny + j;
                    In[q] = I[q] + ry * (I[(size_t)i * Ny + jp] - 2.0 * I[q] + I[(size_t)i * Ny + jm]) + 0.5 * dt * S[q] * F[q];
                    Sn[q] = S[q] - 0.5 * dt * S[q] * F[q];
                }
                tx.solve(&In[j], Ny, tmp.data());
            }
        });

        compute_F(In, Fn);

        // ---- half step 2: explicit in x, implicit in y ----
        tbb::parallel_for(tbb::blocked_range<int>(0, Nx), [&](const tbb::blocked_range<int>& r) {
            std::vector<double> tmp(Ny);
            for (int i = r.begin(); i < r.end(); ++i) {
                int ip = (i == Nx - 1) ? (consistent_bc ? i : i - 1) : i + 1;
                int im = (i == 0) ? (consistent_bc ? i : i + 1) : i - 1;
                for (int j = 0; j < Ny; ++j) {
                    size_t q = (size_t)i * Ny + j;
                    I[q] = In[q] + rx * (In[(size_t)ip * Ny + j] - 2.0 * In[q] + In[(size_t)im * Ny + j]) + 0.5 * dt * Sn[q] * Fn[q];
                    S[q] = Sn[q] - 0.5 * dt * Sn[q] * Fn[q];
                }
                ty.solve(&I[(size_t)i * Ny], 1, tmp.data());
            }
        });

        if (k % out_stride == 0) {
            observe(k);
            if (std::isnan(massI)) { fprintf(stderr, "NaN at t=%g\n", k * dt); return 2; }
            if (k % (out_stride * 200) == 0) printf("  t=%7.1f massI=%.4f massS=%.4f maxI=%.4f area(%g)=%.2f\n", k * dt, massI, massS, maxI, thetas[1 < nth ? 1 : 0], final_area[1 < nth ? 1 : 0]);
        }
    }

    // ---- summary ----
    std::string sum_file = out + "/summary.txt";
    FILE* fsum = fopen(sum_file.c_str(), "w");
    fprintf(fsum, "model %s\nL %g\nD %g\ndx %g\ndt %g\nT %g\nNx %d\nr0 %g\neps %g\nk1 %g\n", model.c_str(), L, D, dx, dt, T, Nx, r0, eps, k1);
    fprintf(fsum, "massI_final %.10f\nmassS_final %.10f\nmaxI_final %.10f\nPR_final %.10f\nvar_final %.12e\n", massI, massS, maxI, PR, var);
    for (size_t a = 0; a < nth; ++a)
        fprintf(fsum, "theta %g peak_area %.6f peak_time %.4f final_area %.6f\n", thetas[a], peak_area[a], peak_time[a], final_area[a]);
    fclose(fsum);
    printf("done: %s\n", out.c_str());
    return 0;
}
