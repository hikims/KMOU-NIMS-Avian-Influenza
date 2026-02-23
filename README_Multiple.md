# Reaction--Diffusion Epidemic Simulation (Multiple Parameter Version)

This repository contains C++ simulation code for a spatial
reaction--diffusion epidemic model under multiple diffusion
coefficients.

All source files are GitHub-ready: - Korean comments have been
translated into English. - No hard-coded absolute paths are used. -
Output directories are configurable by the user.

------------------------------------------------------------------------

## Project Structure

Core simulation files:

-   Multiple_main_github.cpp
-   Multiple_parameters_github.h
-   functions_github.cpp
-   functions_github.hpp
-   initial_github.cpp
-   initial_github.hpp
-   Multiple_exec_github.sh

------------------------------------------------------------------------

## Model Overview

The model simulates epidemic spread with:

-   Spatial diffusion (diffusion coefficient `di`)
-   Susceptible and infected area dynamics
-   Multiple parameter sweeps
-   Region-based output files

Each diffusion coefficient generates a separate output directory.

------------------------------------------------------------------------

## Output Directory Configuration

The simulation does not use hard-coded system paths.

Set your desired project root directory using:

export PROJECT_PATH="/your/project/path"

If not set, the default output location is the current directory (`.`).

Simulation outputs will be written to:

\${PROJECT_PATH}/data/

Each diffusion coefficient produces:

\${PROJECT_PATH}/data/di_x.xxx/

------------------------------------------------------------------------

## Compilation

Example (Linux / macOS):

g++ -O2 Multiple_main_github.cpp functions_github.cpp initial_github.cpp
-o simulation

Or use the provided script:

bash Multiple_exec_github.sh

------------------------------------------------------------------------

## Running the Simulation

./simulation

Output files include:

-   Region-based epidemic area data
-   Time series of infected and susceptible regions
-   Per-diffusion parameter results

------------------------------------------------------------------------

## Dependencies

-   C++17 compatible compiler (g++ / clang++)
-   Linux, macOS, or WSL recommended

------------------------------------------------------------------------

## Notes

-   All paths are user-configurable.
-   No absolute user-specific paths are included.
-   Designed for batch simulation across multiple diffusion
    coefficients.

------------------------------------------------------------------------

## License

Please specify your license here (e.g., MIT, BSD, GPL).
