#!/usr/bin/env python3
"""Supplementary Material for CNSNS-D-26-02308 R2.

Writes into ../supplementary/
  S1_before_after_per_run.csv     originally submitted run versus corrected run, for every (L, D_I) of the main sweep
  S1_before_after_summary.tex     summary table (also used in the response letter)
  per_run_summary_<sweep>.csv     per-run summaries of every sweep of the revised manuscript
  supplementary.tex               source of Supplementary_Material.pdf

Inputs: the originally submitted outputs (ORIG, time series and summary written by the submitted code) and the
revised runs under ../data.
"""
import glob, os, re, sys
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data"); SUP = os.path.join(ROOT, "supplementary"); os.makedirs(SUP, exist_ok=True)
ORIG = os.environ.get("ORIG_DATA", os.path.join(DATA, "original_submission"))
TH = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

def classify(pf, ff):
    if pf < 0.02: return "below detection"
    if ff > 0.5: return "saturation"
    if ff < 0.5 * pf: return "transient"
    return "invasion"
def pfit(x, y): p, b = np.polyfit(np.log(x), np.log(y), 1); return p, np.exp(b)
def offsetfit(L, D):
    f = lambda L, c, L0: c * (L - L0) ** 2
    popt, _ = curve_fit(f, L, D, p0=[1e-5, 5.0]); return popt[0], popt[1], 100 * np.max(np.abs(f(L, *popt) - D) / D)
def boundary(Ls, Ds, ff):
    out = {}
    for i, L in enumerate(Ls):
        f = ff[i]
        if np.isnan(f).any() or f[0] > 0.5 or f[-1] <= 0.5: continue
        j = int(np.argmax(f > 0.5)); out[L] = float(np.exp(np.log(Ds[j - 1]) + (0.5 - f[j - 1]) * (np.log(Ds[j]) - np.log(Ds[j - 1])) / (f[j] - f[j - 1])))
    return out

# ---------------- originally submitted runs ----------------
# summary file of the submitted code: di L Nx Ny dx dy peak_area final_area peak_time final_time massI_final massS_final regime theta
orig = {}
for line in open(os.path.join(ORIG, "fig2", "phase_diagram_summary.txt")):
    if line.startswith("#") or not line.strip(): continue
    t = line.split(); D, L = round(float(t[0]), 6), round(float(t[1]), 3)
    orig[(L, D)] = dict(peak=float(t[6]), final=float(t[7]), MI=float(t[10]), MS=float(t[11]), label_submitted=t[12])
N0 = {}; imax_o = {}
for d in glob.glob(os.path.join(ORIG, "fig4", "L_*", "di_*")):
    f = os.path.join(d, "timeseries.txt")
    if not os.path.exists(f): continue
    L = round(float(re.search(r"L_([\d.]+)", d).group(1)), 3); D = round(float(re.search(r"di_([\d.]+)", d).group(1)), 6)
    a = np.loadtxt(f, comments="#")          # t massI massS areaI centerI maxI
    N0[L] = a[0, 1] + a[0, 2]
    if (L, D) in orig and a[-1, 0] >= 399.9: orig[(L, D)]["Imax"] = a[-1, 5]
print("originally submitted runs:", len(orig))

# ---------------- corrected runs (main sweep) ----------------
new = {}
for d in sorted(glob.glob(os.path.join(DATA, "r1", "fixed", "L_*_D_*"))):
    m = re.search(r"L_([\d.]+)_D_([\d.]+)", d); L, D = round(float(m.group(1)), 3), round(float(m.group(2)), 6)
    a = np.loadtxt(os.path.join(d, "timeseries.txt"), comments="#")   # t massI massS maxI centerI PR var areas...
    new[(L, D)] = dict(peak=a[:, 8].max(), final=a[-1, 8], MI=a[-1, 1], MS=a[-1, 2], Imax=a[-1, 3], N0=a[0, 1] + a[0, 2], dMI=a[-1, 1] - a[0, 1], dMS=a[0, 2] - a[-1, 2],
                       drift=np.max(np.abs(a[:, 1] + a[:, 2] - a[0, 1] - a[0, 2])) / (a[0, 1] + a[0, 2]))
Ls = sorted({k[0] for k in new}); Ds = sorted({k[1] for k in new})
I0mass = 0.1 * np.pi * 0.25 ** 2

rows = []; stat = {"orig": {}, "new": {}}
ffo = np.full((len(Ls), len(Ds)), np.nan); ffn = ffo.copy(); pko = ffo.copy(); pkn = ffo.copy()
for i, L in enumerate(Ls):
    for j, D in enumerate(Ds):
        n = new[(L, D)]; o = orig.get((L, D))
        ffn[i, j] = n["final"] / L ** 2; pkn[i, j] = n["peak"]
        if o:
            ffo[i, j] = o["final"] / L ** 2; pko[i, j] = o["peak"]
            loss_o = 1 - (o["MI"] + o["MS"]) / N0[L]; ratio_o = (N0[L] - I0mass - o["MS"]) / (o["MI"] - I0mass)
        rows.append((L, D,
                     o["final"] if o else np.nan, n["final"], classify(o["peak"] / L**2, o["final"] / L**2) if o else "", classify(n["peak"] / L**2, n["final"] / L**2),
                     o["MI"] if o else np.nan, n["MI"], loss_o if o else np.nan, n["drift"], ratio_o if o else np.nan, n["dMS"] / n["dMI"]))
with open(os.path.join(SUP, "S1_before_after_per_run.csv"), "w") as f:
    f.write("L,D_I,final_area_submitted,final_area_corrected,label_submitted_data,label_corrected,M_I_submitted,M_I_corrected,relative_loss_of_S_plus_I_submitted,max_relative_mass_deviation_corrected,S_removed_over_I_gained_submitted,S_removed_over_I_gained_corrected\n")
    for r in rows: f.write(",".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r) + "\n")
R = np.array([[r[2], r[3], r[6], r[7], r[8], r[9], r[10], r[11]] for r in rows], float)
nsat_o = sum(1 for r in rows if r[4] == "saturation"); nsat_n = sum(1 for r in rows if r[5] == "saturation"); n_o = sum(1 for r in rows if r[4])
dco = boundary(Ls, Ds, ffo); dcn = boundary(Ls, Ds, ffn)
Lo = np.array(sorted(dco)); Do = np.array([dco[L] for L in Lo]); Ln = np.array(sorted(dcn)); Dn = np.array([dcn[L] for L in Ln])
peo = pfit(Lo, Do)[0]; pen = pfit(Ln, Dn)[0]; co, L0o, eo = offsetfit(Lo, Do); cn, L0n, en = offsetfit(Ln, Dn)
c2o = np.exp(np.mean(np.log(Do) - 2 * np.log(Lo))); e2o = 100 * np.max(np.abs(c2o * Lo ** 2 - Do) / Do)
c2n = np.exp(np.mean(np.log(Dn) - 2 * np.log(Ln))); e2n = 100 * np.max(np.abs(c2n * Ln ** 2 - Dn) / Dn)
i18 = Ls.index(18.0); i15 = Ls.index(15.0)
a18o = pfit(np.array(Ds), pko[i18])[0]; a18n = pfit(np.array(Ds), pkn[i18])[0]
ok15 = ~np.isnan(pko[i15]); a15o = pfit(np.array(Ds)[ok15], pko[i15][ok15])[0]; a15n = pfit(np.array(Ds), pkn[i15])[0]
Imo = [o["Imax"] for o in orig.values() if "Imax" in o]; Imn = [n["Imax"] for n in new.values()]
summary = [
    ("Runs available / classified as saturation ($\\theta = 0.01$)", f"{n_o} / {nsat_o}", f"{len(new)} / {nsat_n}"),
    ("Susceptible mass removed per unit infected mass gained", f"{np.nanmin(R[:,6]):.2f}--{np.nanmax(R[:,6]):.2f}", f"{R[:,7].min():.3f}--{R[:,7].max():.3f}"),
    ("Largest relative deviation of $\\sum (S+I)$ from its initial value over a run (a loss in the submitted code)", f"{100*np.nanmin(R[:,4]):.1f}--{100*np.nanmax(R[:,4]):.1f}\\%", f"$\\le {np.ceil(1e4*R[:,5].max() - 1e-9)/100:.2f}\\%$"),
    ("Infected mass $M_I(T)$", f"{np.nanmin(R[:,2]):.1f}--{np.nanmax(R[:,2]):.1f}", f"{R[:,3].min():.1f}--{R[:,3].max():.1f}"),
    ("Maximum density $I_{\\max}(T)$", f"{min(Imo):.2f}--{max(Imo):.2f}", f"{min(Imn):.2f}--{max(Imn):.2f}"),
    ("Boundary $10^3 D_c$ at $L = 11$ / $14$ / $18$", f"{1e3*dco[11.0]:.2f} / {1e3*dco[14.0]:.2f} / {1e3*dco[18.0]:.2f}", f"{1e3*dcn[11.0]:.2f} / {1e3*dcn[14.0]:.2f} / {1e3*dcn[18.0]:.2f}"),
    ("Power-law exponent of $D_c(L)$, $L = 11$--$18$", f"{peo:.1f}", f"{pen:.1f}"),
    ("Offset-quadratic fit: $L_0$ (max. deviation)", f"{L0o:.1f} ({eo:.0f}\\%)", f"{L0n:.1f} ({en:.0f}\\%)"),
    ("Pure $c L^2$ fit: max. deviation (rounded up)", f"{np.ceil(e2o - 1e-9):.0f}\\%", f"{np.ceil(e2n - 1e-9):.0f}\\%"),
    ("Area exponent $\\alpha$, fit at $L = 18$ ($L = 15$)", f"{a18o:.2f} ({a15o:.2f})", f"{a18n:.2f} ({a15n:.2f})"),
]
with open(os.path.join(SUP, "S1_before_after_summary.tex"), "w") as f:
    f.write("\\begin{tabular}{p{7.4cm}cc}\n\\toprule\nQuantity (main sweep, $T = 400$) & Submitted code & Corrected code \\\\\n\\midrule\n")
    for a, b, c in summary: f.write(f"{a} & {b} & {c} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
for a, b, c in summary: print(f"{a:70s} | {b:22s} | {c}")
import json; json.dump([list(x) for x in summary], open(os.path.join(SUP, "S1_before_after_summary.json"), "w"), indent=1)

# ---------------- per-run summaries of all sweeps ----------------
def dump(group_dir, name):
    runs = sorted(glob.glob(os.path.join(group_dir, "*L_*_D_*")))
    if not runs: return 0
    n = 0
    with open(os.path.join(SUP, f"per_run_summary_{name}.csv"), "w") as f:
        f.write("run,model,L,Ly,xc,yc,D_I,dx,dt,T,r0,eps,M_I_final,M_S_final,I_max_final,PR_final,max_rel_mass_deviation," + ",".join(f"peak_area_theta_{t:g},final_area_theta_{t:g}" for t in TH) + "\n")
        for d in runs:
            ts = os.path.join(d, "timeseries.txt")
            if not os.path.exists(os.path.join(d, "summary.txt")): continue
            meta = {}; thetas = None
            for line in open(ts):
                if line.startswith("# model="):
                    for kv in line[2:].split(): k, v = kv.split("="); meta[k] = v
                elif line.startswith("# thetas:"): thetas = [float(x) for x in line.split(":")[1].split()]
                elif not line.startswith("#"): break
            a = np.loadtxt(ts, comments="#"); drift = np.max(np.abs(a[:, 1] + a[:, 2] - a[0, 1] - a[0, 2])) / (a[0, 1] + a[0, 2])
            cols = []
            for t in TH:
                if t in thetas: k = 7 + thetas.index(t); cols += [f"{a[:, k].max():.6g}", f"{a[-1, k]:.6g}"]
                else: cols += ["", ""]
            f.write(",".join([os.path.basename(d), meta["model"], meta["L"], meta.get("Ly", meta["L"]), meta.get("xc", "0.5"), meta.get("yc", "0.5"), meta["D"], meta["dx"], meta["dt"], meta["T"], meta["r0"], meta["eps"],
                              f"{a[-1,1]:.8g}", f"{a[-1,2]:.8g}", f"{a[-1,3]:.6g}", f"{a[-1,5]:.6g}", f"{drift:.3g}"] + cols) + "\n"); n += 1
    return n
counts = {}
for name, path in [("main_sweep", "r1/fixed"), ("coarse_sweep", "r1/fixed_hi_theta"), ("observables_run", "r1/fixed_obs"), ("convergence", "r1/conv"), ("local_model", "r1/local"),
                   ("kernel_model", "r1/kernel"), ("long_runs", "r1/long"), ("kernel_refresh", "r1/kerconv"), ("hotspot_r0_0.25", "r1/hotsize/r0_0.25"), ("hotspot_r0_1.0", "r1/hotsize/r0_1.0"),
                   ("hotspot_eps_0.5", "r1/hotsize/eps_0.5"), ("extended_sweep", "ext"), ("extended_checks", "extchk"), ("reference_2d", "ref2d"),
                   ("main_sweep_conservative_rows", "cons"), ("offcentre_corner", "offc/xc0.25_yc0.25"), ("offcentre_side", "offc/xc0.25_yc0.5"), ("rectangle", "rect"),
                   ("homogenization_r0_0.25", "hoth/r0_0.25"), ("homogenization_r0_1.0", "hoth/r0_1.0"), ("large_diffusivity", "extD"), ("large_diffusivity_check_dt_half", "extDchk/dt_half"), ("large_diffusivity_check_dx_half", "extDchk/dx_half")]:
    if path.startswith("r1/") and not os.path.isdir(os.path.join(DATA, "r1")): path = path[3:]     # first-revision runs written directly below data/
    counts[name] = dump(os.path.join(DATA, path), name)
print("per-run summaries:", counts)
# free-space reference
for tag in ("base", "r0_0.25", "r0_1.0", "eps_0.5", "fine"):
    src = os.path.join(DATA, "free", tag, "free_summary.csv")
    if os.path.exists(src):
        import shutil; shutil.copy(src, os.path.join(SUP, f"free_space_reference_{tag}.csv"))
json.dump(counts, open(os.path.join(SUP, "per_run_counts.json"), "w"))
