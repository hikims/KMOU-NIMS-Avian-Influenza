#!/usr/bin/env python3
"""Free-space (unbounded-domain) reference for the hotspot model, CNSNS-D-26-02308 R2.

With the seed and the hotspot centred at the same point the hotspot model is radially symmetric as long as the
walls are not felt, so the solution in an unbounded plane is obtained from the radial problem

    S_t = -S I chi(r),      I_t = D (I_rr + I_r / r) + S I chi(r),      chi(r) = (1 + tanh((r0 - r)/eps)) / 2,

on 0 < r < R with zero flux at r = R, R so large that the density at R is negligible over the horizon.
Discretisation: conservative finite volumes on annuli (cell centres r_j = (j + 1/2) dr), Crank-Nicolson diffusion,
explicit reaction applied as two half steps around the diffusion step; the discrete mass sum_j (S_j + I_j) V_j is
conserved to rounding. The thresholded area is evaluated from the piecewise-linear profile through the cell centres.

Usage:
    python3 radial_free.py                # full D grid -> ../data/free/free_D_*.npz and ../data/free/free_summary.csv
    python3 radial_free.py --D 0.01       # one value
Options: --T 400 --dr 0.02 --dt 0.001 --r0 0.5 --eps 1.0 --procs 16
"""
import argparse, os, sys, time
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data", "free"))
THETAS = np.array([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])


@njit(cache=True)
def area_above(I, r, dr, theta):
    """area of {I >= theta} for the piecewise-linear radial profile through the cell centres (constant inside r[0])."""
    n = I.shape[0]
    a = 0.0
    if I[0] >= theta:
        a += np.pi * r[0] * r[0]
    for j in range(n - 1):
        f0 = I[j] - theta
        f1 = I[j + 1] - theta
        ra = r[j]
        rb = r[j + 1]
        if f0 >= 0.0 and f1 >= 0.0:
            a += np.pi * (rb * rb - ra * ra)
        elif f0 >= 0.0 and f1 < 0.0:
            rc = ra + dr * f0 / (f0 - f1)
            a += np.pi * (rc * rc - ra * ra)
        elif f0 < 0.0 and f1 >= 0.0:
            rc = ra + dr * f0 / (f0 - f1)
            a += np.pi * (rb * rb - rc * rc)
    return a


@njit(cache=True)
def solve(D, T, dr, dt, r0, eps, k1, init_i, leng, R, out_every, thetas):
    n = int(np.round(R / dr))
    r = (np.arange(n) + 0.5) * dr
    V = 2.0 * np.pi * r * dr
    chi = 0.5 * (1.0 + np.tanh((r0 - r) / eps))
    I = np.zeros(n)
    S = np.ones(n)
    for j in range(n):
        if r[j] <= leng:
            I[j] = init_i
            S[j] = 1.0 - init_i
    # Crank-Nicolson matrices: dI_j/dt = D [ rp (I_{j+1}-I_j) - rm (I_j - I_{j-1}) ] / (r_j dr^2)
    lam = 0.5 * D * dt / (dr * dr)
    lo = np.zeros(n)
    up = np.zeros(n)
    for j in range(n):
        rm = j * dr
        rp = (j + 1) * dr
        lo[j] = lam * rm / r[j]
        up[j] = lam * rp / r[j]
    up[n - 1] = 0.0                       # zero flux at r = R
    diag = 1.0 + lo + up
    # Thomas factors for (diag, -lo, -up)
    cp = np.zeros(n)
    m = np.zeros(n)
    m[0] = diag[0]
    cp[0] = -up[0] / m[0]
    for j in range(1, n):
        m[j] = diag[j] + lo[j] * cp[j - 1]
        cp[j] = -up[j] / m[j]
    nsteps = int(np.round(T / dt))
    stride = max(1, int(np.round(out_every / dt)))
    nout = nsteps // stride + 1
    nth = thetas.shape[0]
    out = np.zeros((nout, 6 + nth))
    rhs = np.zeros(n)
    tmp = np.zeros(n)
    mass0 = np.sum((S + I) * V)

    def_idx = 0
    for k in range(nsteps + 1):
        if k % stride == 0:
            mI = np.sum(I * V)
            mS_def = np.sum((1.0 - S) * V)          # susceptible mass removed (finite in the unbounded plane)
            s2 = np.sum(I * I * V)
            out[def_idx, 0] = k * dt
            out[def_idx, 1] = mI
            out[def_idx, 2] = mS_def
            out[def_idx, 3] = np.max(I)
            out[def_idx, 4] = I[n - 1]              # density at the outer boundary (must stay negligible)
            out[def_idx, 5] = np.sum((S + I) * V) - mass0
            for a in range(nth):
                out[def_idx, 6 + a] = area_above(I, r, dr, thetas[a])
            def_idx += 1
        if k == nsteps:
            break
        # reaction half step
        for j in range(n):
            f = 0.5 * dt * k1 * S[j] * I[j] * chi[j]
            I[j] += f
            S[j] -= f
        # diffusion, Crank-Nicolson
        for j in range(n):
            v = (1.0 - lo[j] - up[j]) * I[j]
            if j > 0:
                v += lo[j] * I[j - 1]
            if j < n - 1:
                v += up[j] * I[j + 1]
            rhs[j] = v
        tmp[0] = rhs[0] / m[0]
        for j in range(1, n):
            tmp[j] = (rhs[j] + lo[j] * tmp[j - 1]) / m[j]
        I[n - 1] = tmp[n - 1]
        for j in range(n - 2, -1, -1):
            I[j] = tmp[j] - cp[j] * I[j + 1]
        # reaction half step
        for j in range(n):
            f = 0.5 * dt * k1 * S[j] * I[j] * chi[j]
            I[j] += f
            S[j] -= f
    return out, r, I, S


def run_one(args):
    D, T, dr, dt, r0, eps, tag = args
    R = max(40.0, 6.0 * np.sqrt(4.0 * D * T) + 10.0)
    t0 = time.time()
    out, r, I, S = solve(D, T, dr, dt, r0, eps, 1.0, 0.1, 0.25, R, 0.2, THETAS)
    os.makedirs(os.path.join(DATA, tag), exist_ok=True)
    np.savez_compressed(os.path.join(DATA, tag, f"free_D_{D:.6f}.npz"), ts=out, thetas=THETAS, r=r[::5], I=I[::5], S=S[::5],
                        D=D, T=T, dr=dr, dt=dt, r0=r0, eps=eps, R=R)
    return D, R, out[-1], time.time() - t0


def d_grid():
    g = set(np.round(10 ** np.linspace(np.log10(5e-4), 0, 67), 7))
    g |= {round(0.001 * k, 6) for k in range(1, 11)}
    g |= {round(0.001 * 2 ** (k / 2), 6) for k in range(20)}
    return sorted(g)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--D", type=float, default=None)
    ap.add_argument("--T", type=float, default=400.0)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--dt", type=float, default=0.001)
    ap.add_argument("--r0", type=float, default=0.5)
    ap.add_argument("--eps", type=float, default=1.0)
    ap.add_argument("--tag", type=str, default="base")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--dmax", type=float, default=None, help="restrict the D grid to D <= dmax")
    a = ap.parse_args()
    Ds = [a.D] if a.D is not None else d_grid()
    if a.dmax is not None: Ds = [D for D in Ds if D <= a.dmax]
    solve(0.01, 1.0, 0.05, 0.01, 0.5, 1.0, 1.0, 0.1, 0.25, 10.0, 0.5, THETAS)      # compile
    jobs = [(D, a.T, a.dr, a.dt, a.r0, a.eps, a.tag) for D in Ds]
    from multiprocessing import Pool
    rows = []
    with Pool(a.procs) as pool:
        for D, R, last, el in pool.imap_unordered(run_one, jobs):
            rows.append((D, R, last)); print(f"D={D:.6f} R={R:.0f} M_I={last[1]:.4f} Imax={last[3]:.4f} I(R)={last[4]:.1e} dM={last[5]:.1e} A(0.01)={last[7]:.3f}  {el:.0f}s", flush=True)
    rows.sort(key=lambda x: x[0])
    with open(os.path.join(DATA, a.tag, "free_summary.csv"), "w") as f:
        f.write("D,R,M_I,M_S_removed,I_max,I_at_R,mass_error," + ",".join(f"A_{th:g}" for th in THETAS) + "\n")
        for D, R, last in rows:
            f.write(f"{D:.6f},{R:.1f}," + ",".join(f"{v:.8g}" for v in last[1:]) + "\n")
    print("FREE DONE", len(rows))
