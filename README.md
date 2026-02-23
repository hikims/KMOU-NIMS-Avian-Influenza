[README.md](https://github.com/user-attachments/files/25475050/README.md)
# Reaction–Diffusion Epidemic Simulation (GitHub Release)

This repository provides **C++ simulation code** for a spatial reaction–diffusion epidemic model and
**figure-generation scripts** (Python/Gnuplot) for producing publication-ready plots.

------------------------------------------------------------------------

## Contents

- [Quick start](#quick-start)
- [Project layout](#project-layout)
- [Single-parameter simulation](#single-parameter-simulation)
- [Multiple-parameter simulation](#multiple-parameter-simulation)
- [Figure generation](#figure-generation)
- [Path configuration](#path-configuration)
- [Dependencies](#dependencies)
- [License](#license)

------------------------------------------------------------------------

## Quick start

### 1) Choose an output root directory

```bash
export PROJECT_PATH="/your/project/path"   # default: "."
```

Outputs will be written under:

- `${PROJECT_PATH}/data/`
- `${PROJECT_PATH}/image/` or `${PROJECT_PATH}/paper_image/` (figures)

### 2) Build (Single or Multiple)

Examples (Linux/macOS/WSL):

```bash
# Single
g++ -O2 Single_main_github.cpp functions_github.cpp initial_github.cpp -o simulation_single

# Multiple
g++ -O2 Multiple_main_github.cpp functions_github.cpp initial_github.cpp -o simulation_multiple
```

Or run the provided scripts (if you prefer):

```bash
bash Single_exec_github.sh
bash Multiple_exec_github.sh
```

### 3) Run

```bash
./simulation_single
# or
./simulation_multiple
```

### 4) Make figures

```bash
pip install numpy pandas matplotlib
export PROJECT_PATH="/your/project/path"
python fig3_github.py
python fig4_github.py
python fig6_github.py
```

For gnuplot figures:

```bash
export DATA_PATH="/your/project/path"
gnuplot fig2_github.gpl
gnuplot fig5_github.gpl
```

------------------------------------------------------------------------

## Project layout

### C++ (shared)

- `functions_github.cpp`, `functions_github.hpp` — core model functions  
- `initial_github.cpp`, `initial_github.hpp` — initial conditions / helpers  
- `*_parameters_github.h` — model & numerical parameters  

### Single

- `Single_main_github.cpp`
- `Single_parameters_github.h`
- `Single_exec_github.sh`

### Multiple

- `Multiple_main_github.cpp`
- `Multiple_parameters_github.h`
- `Multiple_exec_github.sh`

### Figures

- `fig2_github.gpl`, `fig5_github.gpl` (gnuplot)
- `fig3_github.py`, `fig4_github.py`, `fig6_github.py` (python)

------------------------------------------------------------------------

## Single-parameter simulation

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

---

## Multiple-parameter simulation

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

---

## Figure generation

# Figure Generation Scripts (GitHub Version)

This directory contains plotting scripts used to generate figures for
the reaction--diffusion epidemic simulation results.

All files are GitHub-ready:

-   Korean comments have been translated into English.
-   Hard-coded absolute paths have been removed.
-   Output directories are user-configurable.
-   Compatible with Linux, macOS, and WSL.

------------------------------------------------------------------------

## Included Files

-   fig2_github.gpl
-   fig3_github.py
-   fig4_github.py
-   fig5_github.gpl
-   fig6_github.py

------------------------------------------------------------------------

# 1. Gnuplot Scripts

## fig2_github.gpl

Generates time-series plots of infected and susceptible areas.

### Requirements

-   gnuplot

### Usage

export DATA_PATH="/your/project/path" gnuplot fig2_github.gpl

Output is written to:

\${DATA_PATH}/

------------------------------------------------------------------------

## fig5_github.gpl

Generates peak infected area vs diffusion coefficient plots.

### Usage

export DATA_PATH="/your/project/path" gnuplot fig5_github.gpl

------------------------------------------------------------------------

# 2. Python Plotting Scripts

## Requirements

pip install numpy pandas matplotlib

------------------------------------------------------------------------

## fig3_github.py

Generates time-series plots of infected and susceptible regions for
multiple diffusion coefficients.

### Usage

export PROJECT_PATH="/your/project/path" python fig3_github.py

Output directory:

\${PROJECT_PATH}/image/

------------------------------------------------------------------------

## fig4_github.py

Generates heatmaps of cluster area evolution.

### Usage

export PROJECT_PATH="/your/project/path" python fig4_github.py

Output directory:

\${PROJECT_PATH}/paper_image/

------------------------------------------------------------------------

## fig6_github.py

Additional figure generation script (custom visualization).

### Usage

export PROJECT_PATH="/your/project/path" python fig6_github.py

------------------------------------------------------------------------

# Notes

-   All paths are user-configurable via environment variables.
-   If no environment variable is set, the current directory (.) is
    used.
-   Ensure simulation output data exists before running these scripts.

------------------------------------------------------------------------

# License

Specify your license here (e.g., MIT, BSD, GPL).

---

## Path configuration

This repository avoids system-specific absolute paths.

### PROJECT_PATH (C++ + Python)

```bash
export PROJECT_PATH="/your/project/path"
```

If not set, scripts default to the current directory (`.`).

### DATA_PATH (gnuplot)

```bash
export DATA_PATH="/your/project/path"
```

---

## Dependencies

- C++17 compatible compiler (g++/clang++)
- Python 3.9+
- Python packages:

```bash
pip install numpy pandas matplotlib
```

- gnuplot (optional, for `.gpl` scripts)

---

## License

Add your preferred license (e.g., MIT / BSD-3 / GPL-3.0) as `LICENSE`.
