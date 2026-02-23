#include "functions_github.hpp"
#include "initial_github.hpp"
#include "Multiple_parameters_github.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <algorithm>
#include <cstdlib>
#include <string>


// Initial infection center coordinates (initialize_multiple_clusters에서 채워짐)
extern double cluster_x[6];
extern double cluster_y[6];

// Project root for outputs.
// Set via environment variable PROJECT_PATH. Defaults to current directory ".".
static const char* PROJECT_PATH = std::getenv("PROJECT_PATH");
static const char* project_path() { return (PROJECT_PATH && PROJECT_PATH[0]) ? PROJECT_PATH : "."; }


// -------------------- Directory utilities --------------------
static void ensure_dir(const char* path) {
    struct stat st;
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0755) != 0 && errno != EEXIST) { perror("mkdir failed"); exit(1); }
    } else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Path exists but is not a directory: %s\n", path); exit(1);
    }
}

// -------------------- Constants / utilities --------------------
static const double I_LOCK_MIN = 5e-3;  // (박스 산정/면적/지름) 임계값
static const double EPS0       = 0.0;

// ★ 클러스터별 고정 스케일(원하는 값으로 조정)
static const double BOX_SCALE_R[6] = {
    2.0, 4.0, 6.0, 8.0, 10.0, 12.0
};

inline bool in_square(double x, double y, double cx, double cy, double half) {
    return (fabs(x - cx) <= half) && (fabs(y - cy) <= half);
}
// Keep only points strictly inside the box( |x-cx| < half, |y-cy| < half )만 유지,
// Boundary(|x-cx| >= half 또는 |y-cy| >= half)는 모두 0으로.
inline void clamp_box_dirichlet(double **Ir, double cx, double cy, double half)
{
    for (int i=0; i<Nx; ++i) {
        double x = i * dx;
        for (int j=0; j<Ny; ++j) {
            double y = j * dy;
            double dx_ = fabs(x - cx);
            double dy_ = fabs(y - cy);
            if (dx_ >= half || dy_ >= half) {
                Ir[i][j] = 0.0;   // 박스 밖 + Boundary 모두 0
            }
        }
    }
}
inline double estimate_half_from_cluster(double **Ir, int rid) {
    double cx = cluster_x[rid], cy = cluster_y[rid];
    double half_needed = 0.0; bool any = false;
    for (int i=0;i<Nx;++i){ double x=i*dx;
        for (int j=0;j<Ny;++j){
            if (Ir[i][j] < I_LOCK_MIN) continue;
            any = true;
            double h = std::max(fabs(x-cx), fabs(j*dy-cy));
            if (h > half_needed) half_needed = h;
        }
    }
    if (!any) return 0.0;
    double maxHalf = std::max(Lx, Ly);
    return std::min(std::max(half_needed, 0.0), maxHalf);
}
static const double BAND_LO = 0.05, BAND_HI = 0.10;
inline double estimate_band_diameter_from_cluster(double **Ir, int rid) {
    double cx = cluster_x[rid], cy = cluster_y[rid], rmax=-1.0;
    for (int i=0;i<Nx;++i){ double x=i*dx; for (int j=0;j<Ny;++j){
        double v = Ir[i][j];
        if (v>=BAND_LO && v<=BAND_HI){ double r=hypot(x-cx,j*dy-cy); if (r>rmax) rmax=r; }
    }}
    if (rmax < 0.0) for (int i=0;i<Nx;++i){ double x=i*dx; for (int j=0;j<Ny;++j){
        if (Ir[i][j] >= 0.01){ double r=hypot(x-cx,j*dy-cy); if (r>rmax) rmax=r; }
    }}
    if (rmax < 0.0) for (int i=0;i<Nx;++i){ double x=i*dx; for (int j=0;j<Ny;++j){
        if (Ir[i][j] >= I_LOCK_MIN){ double r=hypot(x-cx,j*dy-cy); if (r>rmax) rmax=r; }
    }}
    return (rmax<0.0)? 0.0 : 2.0*rmax;
}
inline double shrink_half_avoid_overlap(int rid, double half,
                                        const bool box_on[6], const double box_half[6]) {
    double cx = cluster_x[rid], cy = cluster_y[rid];
    for (int q=0; q<6; ++q) if (q!=rid && box_on[q] && box_half[q] > 0.0) {
        double cqx = cluster_x[q], cqy = cluster_y[q], hq = box_half[q];
        double dx_inf = fabs(cx - cqx) - hq;
        double dy_inf = fabs(cy - cqy) - hq;
        double gap = std::min(dx_inf, dy_inf);
        if (gap < half) half = std::max(0.0, gap);
    }
    return half;
}

int main(int argc, char **argv)
{
    double **S = mallocDouble(Nx, Ny), **S_new = mallocDouble(Nx, Ny), **S_prev = mallocDouble(Nx, Ny);
    double **F = mallocDouble(Nx, Ny), **F_new = mallocDouble(Nx, Ny);
    double **I_total = mallocDouble(Nx, Ny), **I_mid_total = mallocDouble(Nx, Ny);

    double **I_r[6], **I_r_new[6], **I_r_prev[6];
    for (int r=0;r<6;++r){ I_r[r]=mallocDouble(Nx,Ny); I_r_new[r]=mallocDouble(Nx,Ny); I_r_prev[r]=mallocDouble(Nx,Ny); }

    double *A_ix=mallocSingle(Nx), *B_ix=mallocSingle(Nx), *C_ix=mallocSingle(Nx);
    double *A_iy=mallocSingle(Ny), *B_iy=mallocSingle(Ny), *C_iy=mallocSingle(Ny);
    double **L_ix=mallocDouble(Lp,Nx), **L_iy=mallocDouble(Lp,Ny);
    double *U_ix=mallocSingle(Nx), *U_iy=mallocSingle(Ny);

    double **ix=mallocDouble(Ny,Nx), **x_ix=mallocDouble(Ny,Nx), **z_ix=mallocDouble(Ny,Nx);
    double **iy=mallocDouble(Nx,Ny), **x_iy=mallocDouble(Nx,Ny), **z_iy=mallocDouble(Nx,Ny);
    double *it = mallocSingle((iter / cut) + 1);

    double rix = (di * 0.5 * dt) / (dx * dx);
    double riy = (di * 0.5 * dt) / (dy * dy);

    double **I_tmp = mallocDouble(Nx, Ny);
    initialize_multiple_clusters(S, I_tmp);

    auto nearest_r = [&](int i, int j)->int{
        double x=i*dx, y=j*dy; int best=0; double bestd2=1e300;
        for (int r=0;r<6;++r){ double dx_=x-cluster_x[r], dy_=y-cluster_y[r]; double d2=dx_*dx_+dy_*dy_; if (d2<bestd2){bestd2=d2;best=r;} }
        return best;
    };
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j){
        for (int r=0;r<6;++r){ I_r[r][i][j]=EPS0; I_r_prev[r][i][j]=EPS0; }
        if (I_tmp[i][j] > EPS0) I_r[nearest_r(i,j)][i][j] = I_tmp[i][j];
    }
    freeMallocDouble(I_tmp, Nx);

    for (int i=0;i<Nx;i++){ A_ix[i]=1+2*rix; B_ix[i]=-rix; C_ix[i]=-rix; }
    A_ix[0]=1+rix; B_ix[0]=0.0; C_ix[0]=-rix; A_ix[Nx-1]=1+rix; B_ix[Nx-1]=-rix; C_ix[Nx-1]=0.0;
    for (int j=0;j<Ny;j++){ A_iy[j]=1+2*riy; B_iy[j]=-riy; C_iy[j]=-riy; }
    A_iy[0]=1+riy; B_iy[0]=0.0; C_iy[0]=-riy; A_iy[Ny-1]=1+riy; B_iy[Ny-1]=-riy; C_iy[Ny-1]=0.0;

    LUx(A_ix,B_ix,C_ix,L_ix,U_ix,Nx);
    LUy(A_iy,B_iy,C_iy,L_iy,U_iy,Ny);
    tbb::task_arena arena(30);

    ensure_dir((std::string(project_path()) + "/data").c_str());
    ensure_dir((std::string(project_path()) + "/data/rect").c_str());

    char diam_path[1024]; snprintf(diam_path,sizeof(diam_path),"%s/data/diameter.txt", project_path());
    FILE *fdiam=fopen(diam_path,"a"); if(!fdiam){perror("diameter.txt open fail"); exit(1);}
    fprintf(fdiam,"# k  t  rid  D_band  half  x_min y_min x_max y_max  scale  (I>=%.4g)\n", I_LOCK_MIN); fflush(fdiam);

    char area_file_path[1024]; snprintf(area_file_path,sizeof(area_file_path),"%s/data/area_cluster.txt", project_path());
    FILE *farea = fopen(area_file_path, "a"); if (!farea) { perror("area_cluster.txt open fail"); exit(1); }
    fprintf(farea, "# k  t  A0  A1  A2  A3  A4  A5   (I>=%.4g)\n", I_LOCK_MIN); fflush(farea);

    char rectz_path[1024];

    static const long long trigger_k[] = {10000LL, 50000LL, 100000LL, 500000LL, 1000000LL, 1500000LL};
    const int NUM_TRIG = (int)(sizeof(trigger_k)/sizeof(trigger_k[0]));
    int next_trig_idx = 0;

    struct Box { bool on; double half; } box[6];
    for (int r=0;r<6;++r){ box[r].on=false; box[r].half=-1.0; }

    for (long long k=0; k<=iter; ++k){
        for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j) for (int r=0;r<6;++r) I_r_prev[r][i][j]=I_r[r][i][j];

        for (int j=0;j<Ny;j++){
            S[0][j]=S[1][j]; S[Nx-1][j]=S[Nx-2][j];
            for (int r=0;r<6;++r){ I_r[r][0][j]=I_r[r][1][j]; I_r[r][Nx-1][j]=I_r[r][Nx-2][j]; }
        }
        for (int i=0;i<Nx;i++){
            S[i][0]=S[i][1]; S[i][Ny-1]=S[i][Ny-2];
            for (int r=0;r<6;++r){ I_r[r][i][0]=I_r[r][i][1]; I_r[r][i][Ny-1]=I_r[r][i][Ny-2]; }
        }

        if (k==0){
            for (int r=0;r<6;++r){
                double half_min = estimate_half_from_cluster(I_r[r], r);
                double half = BOX_SCALE_R[r] * half_min;  // ★ 클러스터별 스케일
                double half_cap = std::min(std::min(cluster_x[r], Lx - cluster_x[r]),
                                           std::min(cluster_y[r], Ly - cluster_y[r]));
                half = std::min(half, std::max(0.0, half_cap));
                bool on_tmp[6]={0}; double h_tmp[6]={0};
                half = shrink_half_avoid_overlap(r, half, on_tmp, h_tmp);
                box[r].on = true; box[r].half = half;
            }
        }

        for (int r=0;r<6;++r) if (box[r].on && box[r].half>0.0)
            clamp_box_dirichlet(I_r[r], cluster_x[r], cluster_y[r], box[r].half);

        for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j){
            double s=0.0; for (int r=0;r<6;++r) s+=I_r[r][i][j]; I_total[i][j]=s;
        }

        multiple_compute_F(F, I_total);

        {
            tbb::task_arena arena_local(30);
            arena_local.execute([&]{
                tbb::parallel_for(tbb::blocked_range<int>(0,Ny), [&](const tbb::blocked_range<int>& rj){
                    for (int j=rj.begin(); j<rj.end(); ++j){
                        for (int i=0;i<Nx;++i) S_new[i][j] = S[i][j] + dt * (-S[i][j]*F[i][j]);
                        for (int r=0;r<6;++r){
                            for (int i=0;i<Nx;++i){
                                int jp=(j==Ny-1)? j-1 : j+1, jm=(j==0)? j+1 : j-1;
                                ix[j][i] = dfi(I_r[r][i][j], riy, I_r[r][i][jp], I_r[r][i][j], I_r[r][i][jm], S[i][j], F[i][j]);
                            }
                            solx(L_ix,U_ix,&ix[j][0],&z_ix[j][0],&x_ix[j][0],Nx);
                            for (int i=0;i<Nx;++i) I_r_new[r][i][j] = x_ix[j][i];
                        }
                    }
                });
            });
        }

        for (int r=0;r<6;++r) if (box[r].on && box[r].half>0.0)
            clamp_box_dirichlet(I_r_new[r], cluster_x[r], cluster_y[r], box[r].half);

        for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j){
            double s=0.0; for (int r=0;r<6;++r) s+=I_r_new[r][i][j]; I_mid_total[i][j]=s;
        }
        multiple_compute_F(F_new, I_mid_total);

        {
            tbb::task_arena arena_local(30);
            arena_local.execute([&]{
                tbb::parallel_for(tbb::blocked_range<int>(0,Nx), [&](const tbb::blocked_range<int>& ri){
                    for (int i=ri.begin(); i<ri.end(); ++i){
                        for (int r=0; r<6; ++r) {
                            for (int j=0; j<Ny; ++j) {
                                int ip=(i==Nx-1)? i-1 : i+1, im=(i==0)? i+1 : i-1;
                                iy[i][j] = dfi(I_r_new[r][i][j], rix, I_r_new[r][ip][j], I_r_new[r][i][j], I_r_new[r][im][j],
                                               S_new[i][j], F_new[i][j]);
                            }
                            soly(L_iy,U_iy,&iy[i][0],&z_iy[i][0],&x_iy[i][0],Ny);
                            for (int j=0;j<Ny;++j){
                                I_r[r][i][j] = x_iy[i][j];
                                if (r==0) S[i][j] = S_new[i][j] + dt * (-S_new[i][j]*F_new[i][j]);
                            }
                        }
                    }
                });
            });
        }

        for (int r=0;r<6;++r) if (box[r].on && box[r].half>0.0)
            clamp_box_dirichlet(I_r[r], cluster_x[r], cluster_y[r], box[r].half);

        if (next_trig_idx < NUM_TRIG && k == trigger_k[next_trig_idx]) {
            int pick_r = -1; double best_D = -1.0;
            for (int r=0;r<6;++r){
                // 원하면 이미 on이어도 갱신 허용: 아래 줄을 주석 처리
                // if (box[r].on) continue;
                double D = estimate_band_diameter_from_cluster(I_r[r], r);
                if (D > best_D){ best_D = D; pick_r = r; }
            }
            if (pick_r >= 0 && best_D > 0.0){
                double half_min = estimate_half_from_cluster(I_r[pick_r], pick_r);
                double scale = BOX_SCALE_R[pick_r];                 // ★ 클러스터별 스케일
                double half = scale * half_min;

                double half_cap = std::min(std::min(cluster_x[pick_r], Lx - cluster_x[pick_r]),
                                           std::min(cluster_y[pick_r], Ly - cluster_y[pick_r]));
                half = std::min(half, std::max(0.0, half_cap));

                bool on_tmp[6]; double h_tmp[6];
                for (int q=0;q<6;++q){ on_tmp[q]=box[q].on; h_tmp[q]=box[q].half; }
                half = shrink_half_avoid_overlap(pick_r, half, on_tmp, h_tmp);

                box[pick_r].on = true; box[pick_r].half = half;

                double cx=cluster_x[pick_r], cy=cluster_y[pick_r];
                double x1=std::max(0.0,cx-half), x2=std::min(Lx,cx+half);
                double y1=std::max(0.0,cy-half), y2=std::min(Ly,cy+half);
                fprintf(fdiam, "%lld %.9f %d %.9f %.9f %.9f %.9f %.9f %.9f %.3f\n",
                        (long long)k, k*dt, pick_r, best_D, half, x1,y1,x2,y2, scale);
                fflush(fdiam);
                printf("[k=%lld] pick r=%d  D_band=%.6f  half=%.6f  scale=%.2f  rect=[%.3f %.3f %.3f %.3f]\n",
                       (long long)k, pick_r, best_D, half, scale, x1,y1,x2,y2);
            } else {
                printf("[k=%lld] trigger but no eligible cluster.\n", (long long)k);
            }
            ++next_trig_idx;
        }

        for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j){
            double s=0.0; for (int r=0;r<6;++r) s+=I_r[r][i][j]; I_total[i][j]=s;
        }

        if (k%cut==0){
            printf("Calculation time: %.6lf num_file: %lld\n", dt * k, k / cut);
            int z=(int)(k/cut);
            char fn2[1024]; snprintf(fn2,sizeof(fn2),"%s/data/si%d.txt", project_path(), z);
            FILE* fsi2=fopen(fn2,"w"); if(!fsi2){perror("si(final) open fail"); exit(1);}
            for (int i=0;i<Nx;i++){ fprintf(fsi2,"\n");
                for (int j=0;j<Ny;j++) fprintf(fsi2,"%lf %lf %.6lf %.6lf\n", i*dx, j*dy, S[i][j], I_total[i][j]);
            }
            fclose(fsi2);

            snprintf(rectz_path,sizeof(rectz_path),"%s/data/rect/rect%d.txt", project_path(), z);
            FILE* frz2=fopen(rectz_path,"w");
            if (frz2){
                for (int r=0;r<6;++r) if (box[r].on && box[r].half>0.0){
                    double cx=cluster_x[r], cy=cluster_y[r], h=box[r].half;
                    double x1=std::max(0.0,cx-h), x2=std::min(Lx,cx+h);
                    double y1=std::max(0.0,cy-h), y2=std::min(Ly,cy+h);
                    fprintf(frz2,"%.9f %.9f %.9f %.9f\n", x1,y1,x2,y2);
                }
                fclose(frz2);
            }
        }

        double A[6];
        for (int r=0;r<6;++r){
            double sum=0.0;
            for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
                if (I_r[r][i][j] >= I_LOCK_MIN) sum += dx*dy;
            A[r]=sum;
        }
        fprintf(farea, "%lld %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",
                (long long)k, k*dt, A[0], A[1], A[2], A[3], A[4], A[5]);
        fflush(farea);
    }

    printf("\nSimulation finished.\n");

    if (fdiam) fclose(fdiam);
    if (farea) fclose(farea);

    freeMallocDouble(L_ix, Lp); freeMallocDouble(L_iy, Lp);
    freeMallocSingle(U_ix); freeMallocSingle(U_iy);
    freeMallocSingle(A_ix); freeMallocSingle(B_ix); freeMallocSingle(C_ix);
    freeMallocSingle(A_iy); freeMallocSingle(B_iy); freeMallocSingle(C_iy);

    freeMallocDouble(ix, Ny); freeMallocDouble(x_ix, Ny); freeMallocDouble(z_ix, Ny);
    freeMallocDouble(iy, Nx); freeMallocDouble(x_iy, Nx); freeMallocDouble(z_iy, Ny);

    for (int r=0;r<6;++r){ freeMallocDouble(I_r[r],Nx); freeMallocDouble(I_r_new[r],Nx); freeMallocDouble(I_r_prev[r],Nx); }

    freeMallocDouble(S, Nx); freeMallocDouble(S_new, Nx); freeMallocDouble(S_prev, Nx);
    freeMallocDouble(F, Nx); freeMallocDouble(F_new, Nx);
    freeMallocDouble(I_total, Nx); freeMallocDouble(I_mid_total, Nx);
    freeMallocSingle(it);
    return 0;
}