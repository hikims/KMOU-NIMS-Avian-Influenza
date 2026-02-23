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
