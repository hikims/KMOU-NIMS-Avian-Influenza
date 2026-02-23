
// #define Lx (30.0)
// #define Ly (30.0)
// #define Nx (1500)
// #define Ny (1550)

// #define dx (Lx / Nx)
// #define dy (Ly / Ny)
// #define sigma (0.1)
// #define R0 (0.5)
// #define Lp (2)

// // Reaction–diffusion parameters
// #define k1 (0.2)
// #define k2 (0.001)
// #define N (1.0)
// #define init_i (0.1)
// #define leng (0.25)
// #define di (0.01)

// Time step size (considering stability constraints)
// #define dt (0.0002)

// Number of iterations and output intervals
// #define iter (2000000)
// #define cut (1000)

inline constexpr double k1     = 1.0;    // fix (nondimensionalized)
inline constexpr double sigma  = 1.0;    // fix (nondimensionalized)
// inline constexpr double di     = 0.01; // Control key
inline constexpr double Lx     = 15.0;
inline constexpr double Ly     = 15.0;
inline constexpr int    Nx     = 750;
inline constexpr int    Ny     = 750;
inline constexpr int    Lp     = 2;
inline constexpr double dx     = Lx / Nx;
inline constexpr double dy     = Ly / Ny;
inline constexpr double dt     = 0.0002;
inline constexpr double R0     = 0.5;
inline constexpr double leng   = 0.25;
inline constexpr double init_i = 0.1;
inline constexpr int iter      = 2000000;
inline constexpr int cut       = 1000;