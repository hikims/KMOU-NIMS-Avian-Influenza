# Reaction--Diffusion Epidemic Simulation (Single Parameter Version)

This repository contains C++ simulation code for a spatial
reaction--diffusion epidemic model using a single diffusion coefficient
setting.

All source files are prepared for public release:

-   Korean comments have been translated into English.
-   Hard-coded absolute paths have been removed.
-   Output directories are configurable by the user.

------------------------------------------------------------------------

## Project Structure

Core simulation files:

-   Single_main_github.cpp
-   Single_parameters_github.h
-   functions_github.cpp
-   functions_github.hpp
-   initial_github.cpp
-   initial_github.hpp
-   Single_exec_github.sh

------------------------------------------------------------------------

## Model Overview

The model simulates epidemic spread with:

-   Spatial diffusion (diffusion coefficient `di`)
-   Susceptible and infected area dynamics
-   Region-based epidemic area output

This version runs a single parameter configuration instead of sweeping
multiple diffusion coefficients.

------------------------------------------------------------------------

## Output Directory Configuration

The simulation does not use hard-coded system paths.

Set your desired project root directory using:

export PROJECT_PATH="/your/project/path"

If not set, the default output location is the current directory (`.`).

Simulation outputs will be written to:

\${PROJECT_PATH}/data/

------------------------------------------------------------------------

## Compilation

Example (Linux / macOS):

g++ -O2 Single_main_github.cpp functions_github.cpp initial_github.cpp
-o simulation

Or use the provided script:

bash Single_exec_github.sh

------------------------------------------------------------------------

## Running the Simulation

./simulation

Output files include:

-   Region-based epidemic area data
-   Time series of infected and susceptible regions

------------------------------------------------------------------------

## Dependencies

-   C++17 compatible compiler (g++ / clang++)
-   Linux, macOS, or WSL recommended

------------------------------------------------------------------------

## Notes

-   All paths are user-configurable.
-   No absolute user-specific paths are included.
-   Designed for single-parameter simulation experiments.

------------------------------------------------------------------------

## License

Please specify your license here (e.g., MIT, BSD, GPL).
