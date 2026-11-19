# CNSNS revision: "Coverage and Homogenization Crossovers in Nonlinear Spreading on Bounded Domains: Separating Diffusion, Normalization and Confinement" (title of the first two versions: "Diffusion–Confinement Crossover in Nonlinear Spreading on Bounded Domains"

Code, run scripts, run summaries and time series for the revised versions of
J.-H. Kim and M. Kim, *Communications in Nonlinear Science and Numerical Simulation*, Ms. CNSNS-D-26-02308.

| release    | manuscript version            | contents                                                        |
|------------|-------------------------------|-----------------------------------------------------------------|
| `cnsns-r1` | first revision (October 2026) | corrected simulator, sweeps and analysis of the first revision   |
| `cnsns-r2` | second revision (November 2026) | all of the above, plus the free-space reference, the extended sweep, the additional checks, per-run summaries of every sweep, the time series of every run, and the before/after accounting of the correction |

The simulator in the repository root (`Single_*`, `Multiple_*`) belongs to an earlier paper and is **not**
the one used for the CNSNS manuscript.

## Folder layout

```
code/            simulator, radial free-space solver, run scripts, analysis scripts
summaries/       first revision: summary.txt of every run, regime classification
summaries_r2/    second revision: per-run summary CSV of every sweep, before/after table,
                 extended-sweep tables, free-space reference, numbers quoted in the manuscript
timeseries/      time series of every run (one .tar.xz per sweep) and the outputs of the
                 originally submitted code that enter the before/after comparison
Supplementary_Material.pdf   accounting of the correction, checks, list of files
```

## The correction made in the first revision

The code used for the originally submitted manuscript advanced the susceptible equation by a full time step
in each of the two half-sweeps of the ADI scheme, so that `S` was depleted at twice the stated rate and the
discrete mass `S + I` was not conserved. `code/sim.cpp` corrects this: each half-sweep subtracts from `S`
exactly the increment `(dt/2) S F` that it adds to `I`. `Supplementary_Material.pdf` (Section S1) and
`summaries_r2/S1_before_after_per_run.csv` compare the 90 runs of the main sweep before and after the
correction (ratio of susceptible mass removed to infected mass gained: 2.00 before, 1.00 after).

## Model variants (`--model`)

| option   | force of infection `F`            | description                                         |
|----------|-----------------------------------|-----------------------------------------------------|
| `fixed`  | `k1 * I * chi(x; r0, eps)`        | baseline: transmission hotspot anchored at the seed |
| `local`  | `k1 * I`                          | local mass action (delta-kernel limit)              |
| `kernel` | `k1 * (K * I)`, `K ∝ chi`, ∫K = 1 | nonlocal kernel, reflective (method-of-images) extension, FFT |

Equations: `S_t = -S F`, `I_t = D ΔI + S F`, zero-flux boundary, domain `(0, L) x (0, Ly)` (square by default),
seed `I = init_i` on a disk of radius `leng` at the source position (centre by default), `S = 1 - I`.

## Build and run

```bash
cd code
mkdir -p bin
g++ -std=c++17 -O3 -I/opt/homebrew/include -L/opt/homebrew/lib sim.cpp -ltbb -o bin/sim   # needs oneTBB
bin/sim --model fixed --L 15 --D 0.01 --dx 0.04 --dt 0.001 --T 400 --out OUT \
        --thetas 0.005,0.01,0.02,0.05,0.1,0.2,0.3 --snap 0,100,200,400
python3 radial_free.py --D 0.01          # free-space (unbounded plane) reference; needs numpy, numba
```

Outputs per run: `timeseries.txt` (t, M_I, M_S, I_max, I at the source, participation ratio, spatial variance,
A_I(t; θ) for each θ), `summary.txt` (peak/final areas per θ), optional binary snapshots of `I` and `S`.

Defaults: `r0 = 0.5`, `eps = 1`, `k1 = 1`, `init_i = 0.1`, `leng = 0.25`, output every 0.2 time units,
kernel refreshed every 10 steps (`--conv_every`). Options added for the second revision (the default
behaviour is unchanged, bit for bit): `--Ly <side>` (rectangle, hotspot model), `--xc`, `--yc` (source
position as fractions of the sides).

Boundary rows: by default the implicit sweeps use the `(-1, 1)` zero-flux rows and the explicit sweeps the
mirrored neighbour (`(-2, 2)` rows); this is the production scheme of the main sweep, with a boundary mass
residual below 1e-3. `--consistent_bc` uses `(-1, 1)` rows in both sweeps and conserves the discrete mass to
1e-12; all runs added in the second revision use it, and the main sweep was recomputed with it as a check.

## Scripts

| script                  | runs / output                                                        | used for (second revision) |
|-------------------------|----------------------------------------------------------------------|----------------------------|
| `run_sweep.py`          | hotspot 9×10 grid (dx 0.04, T 400); local & kernel grids (dx 0.02); convergence; long runs | Figs. 2, 6–10, Table 1 |
| `run_extra.py`          | hotspot grid at dx 0.08 with θ up to 0.8; (15, 0.01) with 12 thresholds | Figs. 4, 8 |
| `run_extra2.py`         | more convergence cases; hotspot sizes (r0 = 0.25, 1.0; eps = 0.5); kernel refresh period | Table 1, Fig. 5 |
| `run_bc_check.py`       | boundary-row comparison (six pairs); kernel refresh at T = 10         | Sec. 4.1 |
| `radial_free.py`        | free-space reference: radial finite-volume solver, 93 diffusivities, four hotspots, refined grid | Sec. 4.3, Figs. 2–6 |
| `run_r2_extended.py`    | extended sweep (L = 10…56, D = 0.001…0.72; pairs listed in `ext_needed.json`), two reference runs at L = 40, five production-resolution checks | Sec. 4.3, Fig. 3 |
| `run_r2_robust.py`      | main sweep with mass-conserving rows; off-centre source; rectangles; hotspot sizes for the homogenization crossover | Sec. 4.3 (scope), Supplementary S3 |
| `run_r2_largeD.py`      | L = 44, 48, 56 at diffusivities beyond the extended sweep (D up to 5.8), to locate the homogenization crossover at large D | Sec. 4.3, Supplementary S3 |
| `analyze.py`, `analyze_extra.py` | first-revision analysis: classification, fits, figures, `numbers.tex` | figures unchanged since the first revision |
| `analyze_r2.py`         | free-space comparison, wall effect, extended regime map, homogenization crossover, redrawn Figs. 2–6, `numbers_r2.tex` | Secs. 4.2–4.6 |
| `analyze_r2_robust.py`  | additional checks, `numbers_r2_robust.tex`, Table S2                  | Sec. 4.3 (scope), Supplementary S3 |
| `build_supplementary.py`| before/after table and per-run summary CSV files                      | Supplementary S1, S4 |

The run scripts write to `../data` relative to `code/`; the analysis scripts read from there (the
first-revision sweeps either in `../data` or in `../data/r1`) and write figures and numbers to `../figures`,
`../manuscript` and `../supplementary` (create these folders). To use the archived time series instead of
re-running, unpack `timeseries/*.tar.xz` in the folder that contains `code/`: they restore the `data/` tree.
Total cost of all sweeps: of the order of fifteen hours on a 32-core workstation.

## Classification rules (Eq. (29) of the second revision)

With `a_peak = max_t A_I(t; θ)/L²` and `a_T = A_I(T; θ)/L²`, evaluated in this order:
below detection if `a_peak < 0.02` (called "collapse" in the first revision); saturation if `a_T > 0.5`;
transient if `a_T < a_peak/2`; invasion otherwise. The labels describe the thresholded area only.
The coverage boundary `D_c(L; θ)` is the value of `D` at which `a_T` crosses 0.5, interpolated linearly in
`log D` between the bracketing grid points; the homogenization crossover `D_h(L)` is defined in the same way
from the participation ratio, `PR(T) = 0.9`.
