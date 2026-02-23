#include "functions_github.hpp"
#include "initial_github.hpp"
#include "Single_parameters_github.h"

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
// Project root directory
// Set PROJECT_PATH environment variable to change the output location.
#ifndef PROJECT_PATH
#define PROJECT_PATH "."
#endif
static void ensure_dir(const char* path) {
    struct stat st;
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0755) != 0 && errno != EEXIST) {
            perror("Fail to make directory");
            exit(1);
        }
    } else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Wrong directory path: %s\n", path);
        exit(1);
    }
}

int main(int argc, char **argv)
{
    double **S = mallocDouble(Nx, Ny);
    double **I = mallocDouble(Nx, Ny);
    double **S_new = mallocDouble(Nx, Ny);
    double **I_new = mallocDouble(Nx, Ny);

    double *A_ix = mallocSingle(Nx);
    double *B_ix = mallocSingle(Nx);
    double *C_ix = mallocSingle(Nx);

    double *A_iy = mallocSingle(Ny);
    double *B_iy = mallocSingle(Ny);
    double *C_iy = mallocSingle(Ny);

    double **x_ix = mallocDouble(Ny, Nx);
    double **x_iy = mallocDouble(Nx, Ny);
    double **z_ix = mallocDouble(Ny, Nx);
    double **z_iy = mallocDouble(Nx, Ny);

    double **ix = mallocDouble(Ny, Nx);
    double **iy = mallocDouble(Nx, Ny);

    double **F = mallocDouble(Nx, Ny);
    double **F_new = mallocDouble(Nx, Ny);

    double **L_ix = mallocDouble(Lp, Nx);
    double **L_iy = mallocDouble(Lp, Ny);
    double *U_ix = mallocSingle(Nx);
    double *U_iy = mallocSingle(Ny);
    double *it = mallocSingle((iter / cut) + 1);

    double di = atof(argv[1]);

    double rix = (di * 0.5 * dt) / pow(dx, 2);
    double riy = (di * 0.5 * dt) / pow(dy, 2);

    // -------------------------------------
    // Initial condition
    // -------------------------------------
    // initialize_multiple_clusters(S, I);
    double center_x = Lx / 2;
    double center_y = Ly / 2;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double x = dx * i;
            double y = dy * j;
            I[i][j] = 0.0;
            if (pow(x - center_x, 2) + pow(y - center_y, 2) <= pow(leng, 2))
                I[i][j] = init_i;
            S[i][j] = 1.0 - I[i][j];
        }
    }

    // -------------------------------------
    // LU decomposition setup (with Neumann BC)
    // -------------------------------------
    for (int i = 0; i < Nx; i++) {
        A_ix[i] = 1 + 2 * rix;
        B_ix[i] = -rix;
        C_ix[i] = -rix;
    }
    // Neumann BC in x-direction
    A_ix[0] = 1 + rix;
    B_ix[0] = 0.0;
    C_ix[0] = -rix;

    A_ix[Nx - 1] = 1 + rix;
    B_ix[Nx - 1] = -rix;
    C_ix[Nx - 1] = 0.0;

    for (int j = 0; j < Ny; j++) {
        A_iy[j] = 1 + 2 * riy;
        B_iy[j] = -riy;
        C_iy[j] = -riy;
    }
    // Neumann BC in y-direction
    A_iy[0] = 1 + riy;
    B_iy[0] = 0.0;
    C_iy[0] = -riy;

    A_iy[Ny - 1] = 1 + riy;
    B_iy[Ny - 1] = -riy;
    C_iy[Ny - 1] = 0.0;

    LUx(A_ix, B_ix, C_ix, L_ix, U_ix, Nx);
    LUy(A_iy, B_iy, C_iy, L_iy, U_iy, Ny);

    // Create a TBB task arena to manage thread usage explicitly
    tbb::task_arena arena;

    // -------------------------------------
    // Time iteration
    // -------------------------------------
    for (int k = 0; k <= iter; k++) {

        // -------------------------------
        // Enforce Neumann Boundary
        // -------------------------------
        for (int j = 0; j < Ny; j++) {
            I[0][j]     = I[1][j];
            I[Nx-1][j]  = I[Nx-2][j];
            S[0][j]     = S[1][j];
            S[Nx-1][j]  = S[Nx-2][j];
        }
        for (int i = 0; i < Nx; i++) {
            I[i][0]     = I[i][1];
            I[i][Ny-1]  = I[i][Ny-2];
            S[i][0]     = S[i][1];
            S[i][Ny-1]  = S[i][Ny-2];
        }

        // -------------------------------------
        // Output
        // -------------------------------------
        if (k % cut == 0){
            printf("di : %.3lf Calculation time: %.6lf num_file: %d\n", di, dt * k, k / cut);

            int z = k / cut;
            it[z] = k * dt;

            // --- Prepare output directories: ${PROJECT_PATH}/data and ${PROJECT_PATH}/data/di_%.3f ---
            char dir_base[1024];
            char dir_di[1024];
            snprintf(dir_base, sizeof(dir_base), PROJECT_PATH "data");
            ensure_dir(dir_base);

            // Include di value in the folder name (3 decimal places)
            snprintf(dir_di, sizeof(dir_di), PROJECT_PATH "data/di_%.3f", di);
            ensure_dir(dir_di);

            // --- Prepare output file paths: ${PROJECT_PATH}/data/di_%.3f/si<z> ---
            char filename1[1024];
            snprintf(filename1, sizeof(filename1), "%s/si%d", dir_di, z);

            FILE *fsi = fopen(filename1, "w");
            if (!fsi) {
                perror("Failed to open file");
                exit(1);
            }

            for (int i = 0; i < Nx; i++) {
                fprintf(fsi, "\n");
                for (int j = 0; j < Ny; j++) {
                    fprintf(fsi, "%lf %lf %.6lf %.6lf\n", i * dx, j * dy, S[i][j], I[i][j]);
                }
            }
            fclose(fsi);

            if (I[Nx/2][Ny/2] > 1.0 || isnan(I[Nx/2][Ny/2])) {
                printf("Unstable at step %d (t = %.6lf): I=%.6lf\n", k, k*dt, I[Nx/2][Ny/2]);
                break;
            }
        }

        // -------------------------------------
        // F computation
        // -------------------------------------
        // multiple_compute_F(F, I);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = dx * i;
                double y = dy * j;
                F[i][j] = k1 * I[i][j] * hyptan(dist(center_x - x, center_y - y));
            }
        }

        // --- x-direction ADI ---
        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Ny), [&](const tbb::blocked_range<int> &r) {
                for (int j = r.begin(); j < r.end(); j++) {
                    for (int i = 0; i < Nx; i++) {
                        int jp = (j == Ny - 1) ? j - 1 : j + 1;
                        int jm = (j == 0) ? j + 1 : j - 1;
                        ix[j][i] = dfi(I[i][j], riy, I[i][jp], I[i][j], I[i][jm], S[i][j], F[i][j]);
                    }
                    solx(L_ix, U_ix, &ix[j][0], &z_ix[j][0], &x_ix[j][0], Nx);

                    for (int i = 0; i < Nx; i++) {
                        I_new[i][j] = x_ix[j][i];
                        S_new[i][j] = S[i][j] + dt * (-S[i][j] * F[i][j]);
                    }
                }
            });
        });

        // --- F_new computation ---
        // multiple_compute_F(F_new, I_new);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = dx * i;
                double y = dy * j;
                F_new[i][j] = k1 * I_new[i][j] * hyptan(dist(center_x - x, center_y - y));
            }
        }

        // --- y-direction ADI ---
        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Nx), [&](const tbb::blocked_range<int> &r) {
                for (int i = r.begin(); i < r.end(); i++) {
                    for (int j = 0; j < Ny; j++) {
                        int ip = (i == Nx - 1) ? i - 1 : i + 1;
                        int im = (i == 0) ? i + 1 : i - 1;
                        iy[i][j] = dfi(I_new[i][j], rix, I_new[ip][j], I_new[i][j], I_new[im][j], S_new[i][j], F_new[i][j]);
                    }
                    soly(L_iy, U_iy, &iy[i][0], &z_iy[i][0], &x_iy[i][0], Ny);

                    for (int j = 0; j < Ny; j++) {
                        I[i][j] = x_iy[i][j];
                        S[i][j] = S_new[i][j] + dt * (-S_new[i][j] * F_new[i][j]);
                    }
                }
            });
        });

        // -------------------------------------
        // Check of Infection 
        // -------------------------------------

        // single cluster
        int a = 0;
        int b = 0;
        int c = 0;
        int d = 0;
        double check_i = 0.0;
        double first_term_i = 0.0;
        double second_term_i = 0.0;
        double third_term_i = 0.0;
        double fourth_term_i = 0.0;
        double fifth_term_i = 0.0;

        double check_s = 0.0;
        double first_term_s = 0.0;
        double second_term_s = 0.0;
        double third_term_s = 0.0;
        double fourth_term_s = 0.0;
        double fifth_term_s = 0.0;

        for (int j = 1; j < Ny; j++)
        {
            first_term_i += I[a][j];
            second_term_i += I[b][j];
            first_term_s += S[a][j];
            second_term_s += S[b][j];
        }
        for (int i = 1; i < Nx; i++)
        {
            third_term_i += I[i][c];
            fourth_term_i += I[i][d];
            third_term_s += S[i][c];
            fourth_term_s += S[i][d];
        }

        for (int j = 1; j < Ny; j++)
        {
            for (int i = 1; i < Nx; i++)
            {
                fifth_term_i += I[i][j];
                fifth_term_s += S[i][j];
            }
        }

        check_i = 0.25 * dx * dy                           //
                    * (I[a][c] + I[b][d] + I[a][d] + I[b][c] //
                       + 2 * first_term_i + 2 * second_term_i //
                       + 2 * third_term_i + 2 * fourth_term_i + 4 * fifth_term_i);

        check_s = 0.25 * dx * dy                           //
                    * (S[a][c] + S[b][d] + S[a][d] + S[b][c] //
                       + 2 * first_term_s + 2 * second_term_s //
                       + 2 * third_term_s + 2 * fourth_term_s + 4 * fifth_term_s);

        // -------------------------------------
        // Save per region data
        // -------------------------------------
        char region_dir[1024];
        snprintf(region_dir, sizeof(region_dir), PROJECT_PATH "./region/di_%.3f", di);
        ensure_dir(PROJECT_PATH "./region");
        ensure_dir(region_dir);

        // --- Open per-di output files ---
        char region_file[1024];
        snprintf(region_file, sizeof(region_file), "%s/region.txt", region_dir);
        FILE *fSumI = fopen(region_file, "a");
        if (!fSumI) {
            perror("Failed to open region file");
            exit(1);
        }

        // --- Save outputs ---
        fprintf(fSumI, "%lf %.12lf %.12lf\n", k * dt, check_i, check_s);
        fclose(fSumI);
    }

    // -------------------------------------
    // Time data save
    // -------------------------------------
    FILE *ft = fopen(PROJECT_PATH "./data/time", "w");
    for (int i = 0; i <= iter / cut; i++){
        fprintf(ft, "%.3lf\n", it[i]);
    }
    fclose(ft);

    printf("\nfinish\n");

    // -------------------------------------
    // Free memory
    // -------------------------------------
    freeMallocDouble(L_ix, Lp);
    freeMallocDouble(L_iy, Lp);

    freeMallocDouble(S, Nx);
    freeMallocDouble(S_new, Nx);
    freeMallocDouble(I, Nx);
    freeMallocDouble(I_new, Nx);

    freeMallocDouble(x_ix, Ny);
    freeMallocDouble(x_iy, Nx);
    freeMallocDouble(z_ix, Ny);
    freeMallocDouble(z_iy, Nx);

    freeMallocDouble(F, Nx);
    freeMallocDouble(F_new, Nx);

    freeMallocSingle(U_ix);
    freeMallocSingle(U_iy);
    freeMallocSingle(A_ix);
    freeMallocSingle(B_ix);
    freeMallocSingle(C_ix);
    freeMallocSingle(A_iy);
    freeMallocSingle(B_iy);
    freeMallocSingle(C_iy);
    freeMallocSingle(it);

    return 0;
}
