#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//----------------------------------
// Allocation
//----------------------------------
double** mallocDouble(int row, int col)
{
    double **Mat = (double**) malloc(sizeof(double*) * row);
    if (!Mat) {
        fprintf(stderr, "Memory allocation failed (row pointers)\n");
        exit(1);
    }

    for (int i = 0; i < row; i++)
    {
        Mat[i] = (double*) malloc(sizeof(double) * col);
        if (!Mat[i]) {
            fprintf(stderr, "Memory allocation failed (row %d)\n", i);
            exit(1);
        }

        for (int j = 0; j < col; j++)
            Mat[i][j] = 0.0;
    }
    return Mat;
}

double* mallocSingle(int col)
{
    double *arr = (double*) malloc(sizeof(double) * col);
    if (!arr) {
        fprintf(stderr, "Memory allocation failed (single)\n");
        exit(1);
    }

    for (int j = 0; j < col; j++)
        arr[j] = 0.0;

    return arr;
}

void freeMallocDouble(double **data, int row)
{
    for (int i = 0; i < row; i++)
        free(data[i]);
    free(data);
}

void freeMallocSingle(double *data)
{
    free(data);
}
