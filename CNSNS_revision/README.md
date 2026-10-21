# CNSNS revision: "Diffusion–Confinement Crossover in Nonlinear Spreading on Bounded Domains"

Code, run scripts and run summaries for the revised version (R1, October 2026) of
J.-H. Kim and M. Kim, *Communications in Nonlinear Science and Numerical Simulation*, Ms. CNSNS-D-26-02308.

All figures of the revised manuscript (Figs. 2–9), the numbers quoted in the text (`numbers.tex`) and the
convergence table were produced with the files in this folder. The simulator in the repository root
(`Single_*`, `Multiple_*`) belongs to an earlier paper and is **not** the one used for the revised CNSNS figures.

## Why a new simulator

While re-running the simulations for the revision we found that the time-stepping of the originally submitted
code advanced the susceptible equation by a full time step in each of the two half-sweeps of the ADI scheme,
so that `S` was depleted at twice the stated rate and the discrete mass `S + I` was not conserved.
`code/sim.cpp` corrects this (the susceptible update uses the same half-steps as the ADI sweeps; the discrete
mass is conserved exactly by the reaction step; the residual drift from the boundary quadrature of the diffusion step is below 1e-3), and every figure of the revised manuscript was regenerated with it.

## Model variants (`--model`)

| option   | force of infection `F`           | description                                         |
|----------|----------------------------------|-----------------------------------------------------|
| `fixed`  | `k1 * I * chi(x; r0, eps)`       | baseline: transmission hotspot anchored at the seed |
| `local`  | `k1 * I`                         | local mass action (delta-kernel limit)              |
| `kernel` | `k1 * (K * I)`, `K ∝ chi`, ∫K = 1 | nonlocal kernel, reflective (method-of-images) extension, FFT |

Equations: `S_t = -S F`, `I_t = D ΔI + S F`, zero-flux boundary, square domain `(0, L)^2`,
seed `I = init_i` on a disk of radius `leng` at the centre, `S = 1 - I`.

## Build and run

```bash
cd code
g++ -std=c++17 -O3 -I/opt/homebrew/include -L/opt/homebrew/lib sim.cpp -ltbb -o sim   # needs oneTBB
./sim --model fixed --L 15 --D 0.01 --dx 0.04 --dt 0.001 --T 400 --out OUT \
      --thetas 0.005,0.01,0.02,0.05,0.1,0.2,0.3 --snap 0,100,200,400
```

Outputs per run: `timeseries.txt` (t, M_I, M_S, I_max, I at centre, participation ratio, spatial variance,
A_I(t; θ) for each θ), `summary.txt` (peak/final areas per θ), optional binary snapshots of `I` and `S`.

Defaults: `r0 = 0.5`, `eps = 1`, `k1 = 1`, `init_i = 0.1`, `leng = 0.25`, output every 0.2 time units,
kernel refreshed every 10 steps (`--conv_every`).

Boundary rows: by default the implicit sweeps use the `(-1, 1)` zero-flux rows and the explicit sweeps the mirrored
neighbour (`(-2, 2)` rows); this is the production scheme of the paper, with a boundary mass residual below 1e-3.
`--consistent_bc` uses `(-1, 1)` rows in both sweeps, which conserves the discrete mass to 1e-12; it is used only
for the boundary-treatment check of Section 4.1 (`run_bc_check.py`).

## Reproducing the sweeps

| script             | runs                                                                 | used for            |
|--------------------|----------------------------------------------------------------------|---------------------|
| `run_sweep.py`     | hotspot 9×10 grid (dx 0.04, T 400); local & kernel 5×9 grids (dx 0.02, T 25); convergence (dx 0.02); long runs | Figs. 2, 5, 6, 8, 9, Table 1 |
| `run_extra.py`     | hotspot grid at dx 0.08 with θ up to 0.8; (15, 0.01) with 12 thresholds | Figs. 3, 7          |
| `run_extra2.py`    | more convergence cases; hotspot-size sensitivity (r0 = 0.25, 1.0; eps = 0.5); kernel refresh-period check | Table 1, Fig. 4     |
| `run_bc_check.py`  | boundary-row comparison (mixed rows vs `--consistent_bc`, six pairs at dx 0.08); kernel refresh at T = 10 on both sides of the kernel boundary | Sec. 4.1            |
| `analyze.py`       | classification (Eq. 26 of the paper), fits, Figs. 2, 3, 5–9, `numbers.tex`, `table_conv.tex` | all                 |
| `analyze_extra.py` | convergence table, mass error, hotspot-size and horizon dependence, Fig. 4 | Sec. 4.1, 4.4       |

The scripts write to `../data` relative to `code/`; `../manuscript` must exist for the figure output
(create it or edit `MS` at the top of the analysis scripts). Total cost of all sweeps: about 4 hours on a
32-core workstation.

The front-speed reference for the local model uses the susceptible density ahead of the front, `S_b = 1`
(`D_c = L^2 / (8 pi k1 S_b T^2)` for the half-area criterion).

## Run summaries

`summaries/` contains the `summary.txt` of every run (peak and final thresholded areas, final masses) and the
regime classification at θ = 0.01 (`regimes_fixed_theta0.01.csv`), so that the regime maps and the fitted
boundaries can be checked without re-running the simulations. Time series and snapshots (≈ 2 GB) are
available from the corresponding author on request.

## Classification rules (Eq. 26 of the revised manuscript)

With `a_peak = max_t A_I(t; θ)/L²` and `a_T = A_I(T; θ)/L²`, evaluated in this order:
collapse if `a_peak < 0.02`; saturation if `a_T > 0.5`; transient if `a_T < a_peak/2`; invasion otherwise.
The crossover boundary `D_c(L; θ)` is the value of `D` at which `a_T` crosses 0.5, interpolated linearly in
`log D` between the bracketing grid points.
