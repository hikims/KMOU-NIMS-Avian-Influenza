#include "functions_github.hpp"      
#include "Multiple_parameters_github.h" 
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

double cluster_x[6];
double cluster_y[6];

// --------------------------------------------------------
// LU Decomposition with Neumann Boundary Condition
// --------------------------------------------------------
void LUx(double *a, double *b, double *c, double **L, double *U, int size)
{
    int i;
    L[0][0] = a[0];

    for (i = 1; i < size; i++)  
        L[1][i] = b[i];

    for (i = 0; i < size - 1; i++)  
    {
        U[i] = c[i] / L[0][i];
        L[0][i + 1] = a[i + 1] - L[1][i + 1] * U[i];
    }
}

void solx(double **L, double *U, double *Init, double *z, double *x, int size)
{
    int i;
    z[0] = Init[0] / L[0][0];
    for (i = 1; i < size; i++)
        z[i] = (Init[i] - L[1][i] * z[i - 1]) / L[0][i];

    x[size - 1] = z[size - 1];
    for (i = size - 2; i >= 0; i--)
        x[i] = z[i] - U[i] * x[i + 1];
}

void LUy(double *a, double *b, double *c, double **L, double *U, int size)
{
    int i;
    L[0][0] = a[0];

    for (i = 1; i < size; i++)  
        L[1][i] = b[i];

    for (i = 0; i < size - 1; i++)  
    {
        U[i] = c[i] / L[0][i];
        L[0][i + 1] = a[i + 1] - L[1][i + 1] * U[i];
    }
}

void soly(double **L, double *U, double *Init, double *z, double *x, int size)
{
    int i;
    z[0] = Init[0] / L[0][0];
    for (i = 1; i < size; i++)
        z[i] = (Init[i] - L[1][i] * z[i - 1]) / L[0][i];

    x[size - 1] = z[size - 1];
    for (i = size - 2; i >= 0; i--)
        x[i] = z[i] - U[i] * x[i + 1];
}

// --------------------------------------------------------
// Threading Building Blocks (tbb)
// --------------------------------------------------------
// hypertangent
double hyptan(double b)
{
    return 0.5 * (tanh(b / sigma) + 1);
}

// distance
double dist(double x, double y)
{
    return R0 - sqrt(pow(x, 2) + pow(y, 2));
}

// diffusion term
double dfi(double i, double r, double x, double y, double z, double s, double f)
{
    return i + r * (x - 2.0 * y + z) + dt * 0.5 * s * f;
}

// --------------------------------------------------------
// Initialization of Clusters
// --------------------------------------------------------
void initialize_multiple_clusters(double **S, double **I)
{
    srand(time(NULL));
    const int num_clusters = 6;
    int idx = 0;

    for (int iy = 0; iy < 2; iy++) {
        for (int ix = 0; ix < 3; ix++) {
            double x_start = (Lx / 3.0) * ix;
            double x_end   = (Lx / 3.0) * (ix + 1);
            double y_start = (Ly / 2.0) * iy;
            double y_end   = (Ly / 2.0) * (iy + 1);

            // Base random draw
            double rand_x = x_start + (x_end - x_start) * ((double)rand() / RAND_MAX);
            double rand_y = y_start + (y_end - y_start) * ((double)rand() / RAND_MAX);

            // Clamp coordinates to the [1, 29] range
            if (rand_x < 1.0) rand_x = 1.0;
            if (rand_x > 29.0) rand_x = 29.0;
            if (rand_y < 1.0) rand_y = 1.0;
            if (rand_y > 29.0) rand_y = 29.0;

            cluster_x[idx] = rand_x;
            cluster_y[idx] = rand_y;
            idx++;
        }
    }

    printf(">>> Random infection centers (1~29 range):\n");
    for (int k = 0; k < num_clusters; k++)
        printf("  cluster %d: (%.3f, %.3f)\n", k+1, cluster_x[k], cluster_y[k]);

    for (int i = 0; i < Nx; i++) {
        double x = dx * i;
        for (int j = 0; j < Ny; j++) {
            double y = dy * j;
            double infected = 0.0;
        
            for (int k = 0; k < num_clusters; k++) {
                double dist2 = (x - cluster_x[k]) * (x - cluster_x[k]) +
                               (y - cluster_y[k]) * (y - cluster_y[k]);
                if (dist2 <= leng * leng) {
                    infected = init_i;
                    break;
                }
            }
            I[i][j] = infected;
            S[i][j] = 1.0 - infected;
        }
    }
}



void initialize_single_cluster(double **S, double **I)
{
    // No randomness needed
    const int num_clusters = 1;

    // Fix the center coordinates
    cluster_x[0] = Lx * 0.5;
    cluster_y[0] = Ly * 0.5;

    printf(">>> Single infection center at domain center:\n");
    printf("  cluster 1: (%.3f, %.3f)\n", cluster_x[0], cluster_y[0]);

    for (int i = 0; i < Nx; i++) {
        double x = dx * i;
        for (int j = 0; j < Ny; j++) {
            double y = dy * j;

            // Squared distance to the center cluster
            double dist2 =
                (x - cluster_x[0]) * (x - cluster_x[0]) +
                (y - cluster_y[0]) * (y - cluster_y[0]);

            double infected = (dist2 <= leng * leng) ? init_i : 0.0;

            I[i][j] = infected;
            S[i][j] = 1.0 - infected;
        }
    }
}

// -----------------------------
// 1) DSU (Union-Find)
// -----------------------------
namespace {
    struct DSU {
        std::vector<int> p, r;
        explicit DSU(int n) : p(n), r(n, 0) { for (int i = 0; i < n; ++i) p[i] = i; }
        int find(int a) { return p[a] == a ? a : p[a] = find(p[a]); }
        void unite(int a, int b) {
            a = find(a); b = find(b); if (a == b) return;
            if (r[a] < r[b]) std::swap(a, b);
            p[b] = a; if (r[a] == r[b]) r[a]++;
        }
    };
    inline double dist2_xy(double x1, double y1, double x2, double y2) {
        double ddx = x1 - x2;
        double ddy = y1 - y2;
        return ddx*ddx + ddy*ddy;
    }
    } // anonymous namespace
    
    namespace {
    MergedClusters merge_close_clusters(const double *cx_in,
                                        const double *cy_in,
                                        int Nc,
                                        double r2)
    {
        DSU dsu(Nc);
        for (int i = 0; i < Nc; ++i) {
            for (int j = i + 1; j < Nc; ++j) {
                if (dist2_xy(cx_in[i], cy_in[i], cx_in[j], cy_in[j]) <= r2) {
                    dsu.unite(i, j);
                }
            }
        }
    
        // root → indices
        std::unordered_map<int, std::vector<int>> groups;
        groups.reserve(Nc);
        for (int i = 0; i < Nc; ++i) groups[dsu.find(i)].push_back(i);
    
        // 센트로이드 계산
        MergedClusters out;
        out.cx.reserve(groups.size());
        out.cy.reserve(groups.size());
        out.size.reserve(groups.size());
    
        for (auto &kv : groups) {
            const auto &v = kv.second;
            double sx = 0.0, sy = 0.0;
            for (int idx : v) { sx += cx_in[idx]; sy += cy_in[idx]; }
            out.cx.push_back(sx / v.size());
            out.cy.push_back(sy / v.size());
            out.size.push_back((int)v.size());
        }
        return out;
    }
    
    inline double merge_radius2() {
        const double r = 3.0 * std::max(dx, dy);
        return r * r;
    }
    } 
    
    void multiple_compute_F(double **F_out, double **I_in)
    {
        const int    Nc = 6;                 // 원래 클러스터 개수(프로젝트 설정에 맞게)
        const double r2 = merge_radius2();   // 병합 반경^2
        MergedClusters mc = merge_close_clusters(cluster_x, cluster_y, Nc, r2);
        const int    Nm = (int)mc.cx.size(); // 병합 후 개수
    
        // 격자 전 범위 병렬
        tbb::parallel_for(
            tbb::blocked_range2d<int>(0, Nx, 0, Ny),
            [&](const tbb::blocked_range2d<int> &br){
                for (int i = br.rows().begin(); i < br.rows().end(); ++i) {
                    const double x = dx * i;
                    for (int j = br.cols().begin(); j < br.cols().end(); ++j) {
                        const double y = dy * j;
    
                        double f_max = 0.0;
                        for (int k = 0; k < Nm; ++k) {
                            // dist()는 프로젝트의 기존 함수: R0 - sqrt(x^2 + y^2)
                            // hyptan()도 기존 함수: 0.5*(tanh(b/sigma)+1)
                            double val = hyptan( dist(mc.cx[k] - x, mc.cy[k] - y) );
                            if (val > f_max) f_max = val;
                        }
                        F_out[i][j] = k1 * I_in[i][j] * f_max;
                    }
                }
            }
        );
    }

void single_compute_F(double **F_out, double **I_in)
{
    // 단일 클러스터 좌표 (이미 어딘가에서 설정되어 있다고 가정)
    const double cx = cluster_x[0];
    const double cy = cluster_y[0];

    tbb::parallel_for(
        tbb::blocked_range2d<int>(0, Nx, 0, Ny),
        [&](const tbb::blocked_range2d<int> &br){
            for (int i = br.rows().begin(); i < br.rows().end(); ++i) {
                const double x = dx * i;
                for (int j = br.cols().begin(); j < br.cols().end(); ++j) {
                    const double y = dy * j;

                    // 단일 클러스터 기여만 계산
                    const double val = hyptan( dist(cx - x, cy - y) );
                    F_out[i][j] = k1 * I_in[i][j] * val;
                }
            }
        }
    );
}