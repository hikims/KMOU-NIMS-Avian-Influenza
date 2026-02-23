
#include <vector>

//---------------------------
// FUNCTION PROTPTYPE
//---------------------------
struct MergedClusters {
    std::vector<double> cx;   
    std::vector<double> cy;   
    std::vector<int>    size; 
};

void LUx(double *a, double *b, double *c, double **L, double *U, int size);
void solx(double **L, double *U, double *Init, double *z, double *x, int size);
void LUy(double *a, double *b, double *c, double **L, double *U, int size);
void soly(double **L, double *U, double *Init, double *z, double *x, int size);
double hyptan(double b);
double dist(double x, double y);
double dfi(double i, double r, double x, double y, double z, double s, double f);
void initialize_multiple_clusters(double **S, double **I);
void initialize_single_cluster(double **S, double **I);
extern double cluster_x[6];
extern double cluster_y[6];

void multiple_compute_F(double **F_out, double **I_in);
void single_compute_F(double **F_out, double **I_in);