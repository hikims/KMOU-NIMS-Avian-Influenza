#!/usr/bin/env python3
"""Second-revision analysis (CNSNS-D-26-02308 R2).

Inputs  ../data/r1/...        first-revision runs (production sweep `fixed`, coarse sweeps, variants, checks)
        ../data/free/<tag>/   free-space (unbounded plane) reference, radial solver radial_free.py
        ../data/ext/          extended (L, D) sweep, coarse resolution, mass-conserving boundary rows
        ../data/ref2d/        two-dimensional runs in a large domain (validation of the radial solver)
Outputs ../figures/*.png|eps (also copied to ../manuscript/), ../manuscript/numbers_r2.tex,
        ../supplementary/*.csv|tex
"""
import glob, os, re, sys, json, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data"); R1 = os.path.join(DATA, "r1") if os.path.isdir(os.path.join(DATA, "r1")) else DATA   # first-revision runs
FIG = os.path.join(ROOT, "figures"); MS = os.path.join(ROOT, "manuscript"); SUP = os.path.join(ROOT, "supplementary")
for d in (FIG, SUP): os.makedirs(d, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 150, "savefig.dpi": 300})
TH = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
T0 = 400.0
numbers = {}
def N(k, v): numbers[k] = v
def up(v, nd):
    """upper bound: |v| rounded UP to nd decimals (used for every number quoted as 'at most', 'within', 'up to')"""
    return f"{np.ceil(abs(v) * 10 ** nd - 1e-9) / 10 ** nd:.{nd}f}"
def fmt_e(x, d=1):
    m, e = f"{x:.{d}e}".split("e"); return rf"{m}\times 10^{{{int(e)}}}"
def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".eps")); fig.savefig(os.path.join(FIG, name + ".png")); plt.close(fig)

# ----------------------------------------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------------------------------------
class Run:
    def __init__(self, d):
        self.dir = d; self.thetas = None; self.meta = {}
        with open(os.path.join(d, "timeseries.txt")) as f:
            for line in f:
                if line.startswith("# thetas:"): self.thetas = np.array([float(x) for x in line.split(":")[1].split()])
                elif line.startswith("# model="):
                    for kv in line[2:].split(): k, v = kv.split("="); self.meta[k] = v
                elif not line.startswith("#"): break
        self.ts = np.loadtxt(os.path.join(d, "timeseries.txt"), comments="#")
        self.L = float(self.meta["L"]); self.D = float(self.meta["D"]); self.model = self.meta["model"]
        self.Nx = int(self.meta["Nx"]); self.dx = float(self.meta["dx"])
    def t(self): return self.ts[:, 0]
    def col(self, th): return 7 + int(np.argmin(np.abs(self.thetas - th)))
    def area(self, th): return self.ts[:, self.col(th)]
    def at(self, T): return int(np.argmin(np.abs(self.t() - T)))
    def fractions(self, th, T):
        a = self.area(th) / self.L ** 2; j = self.at(T); return a[: j + 1].max(), a[j]
    def complete(self, T=T0): return abs(self.t()[-1] - T) < 1e-6 or self.t()[-1] >= T

def load_group(path, need_T=None):
    runs = {}
    for d in sorted(glob.glob(os.path.join(path, "*L_*_D_*"))):
        if not os.path.exists(os.path.join(d, "summary.txt")): continue
        r = Run(d)
        if need_T is not None and r.t()[-1] < need_T - 1e-6: continue
        runs[(round(r.L, 3), round(r.D, 6))] = r
    return runs

class Free:
    """free-space coverage function A_free(T; D, theta) and mass M_free(T; D) from the radial solver"""
    def __init__(self, tag):
        self.ts = {}; self.prof = {}
        for f in glob.glob(os.path.join(DATA, "free", tag, "free_D_*.npz")):
            z = np.load(f); D = round(float(z["D"]), 6); self.ts[D] = z["ts"]; self.prof[D] = (z["r"], z["I"])
        self.Dg = np.array(sorted(self.ts))
    def _row(self, D, T): ts = self.ts[D]; return ts[int(np.argmin(np.abs(ts[:, 0] - T)))]
    def A_grid(self, T=T0, th=0.01): return np.array([self._row(D, T)[6 + TH.index(th)] for D in self.Dg])
    def M_grid(self, T=T0): return np.array([self._row(D, T)[1] for D in self.Dg])
    def A(self, D, T=T0, th=0.01): return np.interp(np.log(D), np.log(self.Dg), self.A_grid(T, th))
    def M(self, D, T=T0): return np.interp(np.log(D), np.log(self.Dg), self.M_grid(T))
    def Dc(self, L, T=T0, th=0.01, level=0.5, grid=None):
        """smallest D at which A_free = level * L^2 (log-linear interpolation on `grid`, default the fine D grid)"""
        Dg = self.Dg if grid is None else np.asarray(grid); A = self.A_grid(T, th) if grid is None else np.array([self.A(D, T, th) for D in Dg])
        tgt = level * L * L
        if A[0] > tgt: return np.nan
        idx = np.where(A > tgt)[0]
        if len(idx) == 0: return np.nan
        j = idx[0]
        return float(np.exp(np.log(Dg[j - 1]) + (tgt - A[j - 1]) * (np.log(Dg[j]) - np.log(Dg[j - 1])) / (A[j] - A[j - 1])))

def disk_in_square(A, L):
    """area of the centred disk of area A that lies inside the square of side L"""
    r = np.sqrt(max(A, 0.0) / np.pi); h = L / 2
    if r <= h: return A
    if r >= h * np.sqrt(2): return L * L
    seg = r * r * np.arccos(h / r) - h * np.sqrt(r * r - h * h)
    return np.pi * r * r - 4 * seg

def classify(peak_frac, final_frac):
    if peak_frac < 0.02: return "below detection"
    if final_frac > 0.50: return "saturation"
    if final_frac < 0.5 * peak_frac: return "transient"
    return "invasion"
REGIMES = ["below detection", "transient", "invasion", "saturation"]
RCOL = {"below detection": "#d9d9d9", "transient": "#fdb462", "invasion": "#80b1d3", "saturation": "#8dd3c7"}

def regime_table(runs, th, T):
    Ls = sorted({k[0] for k in runs}); Ds = sorted({k[1] for k in runs})
    lab = np.full((len(Ls), len(Ds)), "", dtype=object); ff = np.full((len(Ls), len(Ds)), np.nan); pk = ff.copy()
    for (L, D), r in runs.items():
        pf, fr = r.fractions(th, T); i, j = Ls.index(L), Ds.index(D); lab[i, j] = classify(pf, fr); ff[i, j] = fr; pk[i, j] = pf * L * L
    return np.array(Ls), np.array(Ds), lab, ff, pk

def crossing(Ds, f, level=0.5):
    """first upward crossing of `level` (log-linear in D); returns (D_c, flag) with flag -1/0/+1 = below grid/inside/none"""
    ok = ~np.isnan(f); d, f = np.asarray(Ds)[ok], np.asarray(f)[ok]
    if len(d) == 0: return (np.nan, 1)
    if f[0] > level: return (d[0], -1)
    idx = np.where(f > level)[0]
    if len(idx) == 0: return (d[-1], +1)
    j = idx[0]
    return (float(np.exp(np.log(d[j - 1]) + (level - f[j - 1]) * (np.log(d[j]) - np.log(d[j - 1])) / (f[j] - f[j - 1]))), 0)

def powerfit(x, y):
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float)); p, b = np.polyfit(lx, ly, 1)
    return p, np.exp(b)
def offsetfit(L, D):
    f = lambda L, c, L0: c * (L - L0) ** 2
    popt, _ = curve_fit(f, L, D, p0=[1e-5, 5.0]); return popt[0], popt[1], 100 * np.max(np.abs(f(L, *popt) - D) / D)

def draw_map(ax, Ls, Ds, lab, title=None):
    Z = np.array([[REGIMES.index(l) if l else np.nan for l in row] for row in lab], float)
    lx = np.log10(Ds); e = np.zeros(len(Ds) + 1); e[1:-1] = 0.5 * (lx[:-1] + lx[1:]); e[0] = lx[0] - (e[1] - lx[0]); e[-1] = lx[-1] + (lx[-1] - e[-2])
    ye = np.zeros(len(Ls) + 1); ye[1:-1] = 0.5 * (Ls[:-1] + Ls[1:]); ye[0] = Ls[0] - (ye[1] - Ls[0]); ye[-1] = Ls[-1] + (Ls[-1] - ye[-2])
    ax.pcolormesh(10 ** e, ye, Z, cmap=ListedColormap([RCOL[r] for r in REGIMES]), norm=BoundaryNorm(np.arange(-0.5, 4.5, 1), 4), shading="flat", edgecolors="w", linewidth=0.3)
    ax.set_xscale("log"); ax.set_xlabel("Diffusion coefficient $D_I$"); ax.set_ylabel("Domain size $L$"); ax.set_yticks(Ls)
    if title: ax.set_title(title, fontsize=9)

# ==========================================================================================================
free = Free("base")
fixed = load_group(os.path.join(R1, "fixed"))
print(f"free-space D values: {len(free.Dg)}; production runs: {len(fixed)}")
Ls, Ds, lab, ff, pk = regime_table(fixed, 0.01, T0)

# ---------- 1. validation of the radial solver against 2-D runs in which the walls are not felt ----------
val = []
for (L, D), r in fixed.items():
    rr, II = free.prof[D] if D in free.prof else (None, None)
    if rr is None: continue
    if np.interp(L / 2, rr, II) < 1e-8:                      # free-space density at the nearest wall negligible
        val.append((L, D, r.area(0.01)[-1], free.A(D), r.ts[-1, 1], free.M(D)))
val = np.array(val)
dA = 100 * np.abs(val[:, 2] / val[:, 3] - 1); dM = 100 * np.abs(val[:, 4] / val[:, 5] - 1)
print(f"radial vs 2-D production, wall-free runs: n={len(val)}, max |dA|={dA.max():.3f}%, max |dM|={dM.max():.3f}%")
N("valN", len(val)); N("valA", up(dA.max(), 2)); N("valM", up(dM.max(), 2))
ref2d = load_group(os.path.join(DATA, "ref2d"), need_T=T0)
for (L, D), r in sorted(ref2d.items()):
    print(f"   ref2d L={L:g} D={D:g}: A(0.01) 2-D {r.area(0.01)[-1]:.3f} radial {free.A(D):.3f} ({100*(r.area(0.01)[-1]/free.A(D)-1):+.3f}%), M_I {r.ts[-1,1]:.4f} vs {free.M(D):.4f}")
if not ref2d: N("refA", "--"); N("refM", "--")
if ref2d:
    N("refA", up(max(100*abs(r.area(0.01)[-1]/free.A(D)-1) for (L,D),r in ref2d.items()), 2))
    N("refM", up(max(100*abs(r.ts[-1,1]/free.M(D)-1) for (L,D),r in ref2d.items()), 2))

# ---------- 2. production boundary versus the free-space level curve ----------
dc_conf = {L: crossing(Ds, ff[i]) for i, L in enumerate(Ls)}
Lf = np.array([L for L in Ls if dc_conf[L][1] == 0]); Df = np.array([dc_conf[L][0] for L in Lf])
dc_free_same = np.array([free.Dc(L, grid=Ds) for L in Lf])        # same interpolation grid as the confined sweep
dc_free_fine = np.array([free.Dc(L) for L in Lf])
rel = 100 * (Df / dc_free_same - 1)
print("L, D_c confined, D_c free (same grid), rel %:", [(int(L), round(a, 5), round(b, 5), round(c, 2)) for L, a, b, c in zip(Lf, Df, dc_free_same, rel)])
N("bdyRel", up(np.max(np.abs(rel)), 1))
# wall effect on the coverage at the grid points bracketing the boundary and over the whole sweep
delta = np.full(ff.shape, np.nan)
for i, L in enumerate(Ls):
    for j, D in enumerate(Ds): delta[i, j] = ff[i, j] - disk_in_square(free.A(D), L) / L ** 2
near = []
for i, L in enumerate(Ls):
    if dc_conf[L][1] != 0: continue
    j = int(np.argmax(ff[i] > 0.5)); near += [abs(delta[i, j - 1]), abs(delta[i, j])]
print(f"wall effect at the points bracketing the boundary: max {100*max(near):.2f}% of the domain; over the sweep: max {100*np.nanmax(delta):.2f}%, min {100*np.nanmin(delta):.2f}%")
N("wallNear", up(100*max(near), 1)); N("wallMax", up(100*np.nanmax(delta), 1))
im, jm = np.unravel_index(np.nanargmax(delta), delta.shape); N("wallMaxL", f"{Ls[im]:g}"); N("wallMaxD", f"{Ds[jm]:g}"); N("wallMaxCov", f"{ff[im, jm]:.2f}")
# coverage at which the wall effect first exceeds one per cent of the domain
cov1 = [ff[i, j] for i in range(len(Ls)) for j in range(len(Ds)) if delta[i, j] > 0.01]
N("wallOnsetCov", f"{min(cov1):.2f}" if cov1 else "--")
dMass = np.array([[100 * (fixed[(L, D)].ts[-1, 1] / free.M(D) - 1) for D in Ds] for L in Ls]); N("wallMass", up(np.max(np.abs(dMass)), 1))

# ---------- 3. free-space coverage function: maximum, limiting domain size, local slope ----------
Ag = free.A_grid(); Mg = free.M_grid(); imax = int(np.argmax(Ag)); Lstar = np.sqrt(2 * Ag[imax])
N("AfreeMax", f"{Ag[imax]:.0f}"); N("DfreeMax", f"{free.Dg[imax]:.2f}"); N("Lstar", f"{Lstar:.1f}"); N("MfreeAtMax", f"{Mg[imax]:.0f}")
N("LstarHeur", f"{np.sqrt(2 * Mg[imax] / (np.e * 0.01)):.0f}")
print(f"free space: max A={Ag[imax]:.1f} at D={free.Dg[imax]:.4f}, L*={Lstar:.2f}; heuristic sqrt(2M/(e theta))={np.sqrt(2*Mg[imax]/(np.e*0.01)):.1f}; M_I there {Mg[imax]:.1f}")
sl = np.gradient(np.log(Ag), np.log(free.Dg))
m1 = (free.Dg >= 0.001) & (free.Dg <= 0.01)
N("slopeMin", f"{sl[m1].min():.2f}"); N("slopeMax", f"{sl[m1].max():.2f}"); print(f"local slope of A_free over 1e-3..1e-2: {sl[m1].min():.2f}..{sl[m1].max():.2f}")
N("slopeAtFive", f"{np.interp(np.log(0.05), np.log(free.Dg), sl):.2f}")
for th in (0.005, 0.02, 0.05):
    A_ = free.A_grid(T0, th); N("Lstar" + {0.005: "Low", 0.02: "Two", 0.05: "Five"}[th], f"{np.sqrt(2 * A_.max()):.1f}")
for T in (200.0, 300.0):
    A_ = free.A_grid(T, 0.01); N("Lstar" + ("TwoHundred" if T == 200 else "ThreeHundred"), f"{np.sqrt(2 * A_.max()):.1f}")
# R1 fits and their extrapolation against the free-space curve
c_off, L0, err_off = offsetfit(Lf, Df); pexp, pref = powerfit(Lf, Df)
ext_err = {}
for L in (20, 24, 28, 32, 36):
    dcf = free.Dc(L); ext_err[L] = 100 * (c_off * (L - L0) ** 2 / dcf - 1) if not np.isnan(dcf) else np.nan
    print(f"   L={L}: free D_c={dcf:.5f}, offset-quadratic extrapolation {c_off*(L-L0)**2:.5f} ({ext_err[L]:+.0f}%), power law {pref*L**pexp:.5f} ({100*(pref*L**pexp/dcf-1):+.0f}%)")
N("extErrTwentyFour", f"{abs(ext_err[24]):.0f}"); N("extErrTwentyEight", f"{abs(ext_err[28]):.0f}"); N("extErrThirtyTwo", f"{abs(ext_err[32]):.0f}"); N("extErrThirtySix", f"{abs(ext_err[36]):.0f}")
N("DcFreeTwentyFour", f"{free.Dc(24):.4f}"); N("DcFreeThirtyTwo", f"{free.Dc(32):.4f}"); N("DcFreeThirtySix", f"{free.Dc(36):.3f}")
N("powErrThirtyTwo", f"{abs(100*(pref*32**pexp/free.Dc(32)-1)):.0f}")

# ---------- Fig. 2: regime map with the free-space level curve ----------
fig, ax = plt.subplots(figsize=(6.4, 4.4)); draw_map(ax, Ls, Ds, lab)
ax.plot(Df, Lf, "ko", ms=4, label="numerical boundary $D_c(L)$")
Lc2 = [L for L in Ls if dc_conf[L][1] != 0]
if Lc2: ax.plot([dc_conf[L][0] for L in Lc2], Lc2, "o", mfc="w", mec="k", ms=4, label="boundary outside grid")
Lline = np.linspace(Ls[0] - 0.5, Ls[-1] + 0.5, 80)
ax.plot([free.Dc(L) for L in Lline], Lline, "-", color="C3", lw=1.8, label=r"free-space reference $A_{\rm free}(T;D_I) = L^2/2$")
ax.plot(c_off * (Lline - L0) ** 2, Lline, "k-", lw=1.0, label=rf"$D_c = c\,(L - L_0)^2$, $L_0 = {L0:.1f}$")
ax.plot(pref * Lline ** pexp, Lline, "k--", lw=1.0, label=rf"$D_c \propto L^{{{pexp:.1f}}}$")
ax.set_xlim(10 ** (np.log10(Ds[0]) - 0.06), 10 ** (np.log10(Ds[-1]) + 0.06)); ax.set_ylim(Ls[0] - 0.5, Ls[-1] + 0.5)
h, l = ax.get_legend_handles_labels(); h += [Patch(color=RCOL["invasion"], label="invasion"), Patch(color=RCOL["saturation"], label="saturation")]
ax.legend(handles=h, loc="upper left", fontsize=7); ax.set_title(r"Regime map, hotspot model ($\theta = 0.01$, $T = 400$)", fontsize=9)
fig.tight_layout(); save(fig, "regime_map")

# ---------- Fig. (threshold sensitivity): same content as before, label renamed ----------
hi = load_group(os.path.join(R1, "fixed_hi_theta"))
TH_PROD = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]
fig = plt.figure(figsize=(7.2, 7.0)); gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.95], hspace=0.45, left=0.1, right=0.97, top=0.9, bottom=0.07)
gsa = gs[0].subgridspec(1, 4, wspace=0.12)
for k, (th, grp) in enumerate([(0.005, fixed), (0.05, fixed), (0.3, hi), (0.6, hi)]):
    ax = fig.add_subplot(gsa[0, k]); L_, D_, lab_, ff_, _ = regime_table(grp, th, T0); draw_map(ax, L_, D_, lab_, title=rf"$\theta = {th}$")
    dd = {L: crossing(D_, ff_[i]) for i, L in enumerate(L_)}
    Lc = [L for L in L_ if dd[L][1] == 0]; ax.plot([dd[L][0] for L in Lc], Lc, "ko", ms=4)
    Lo = [L for L in L_ if dd[L][1] != 0]; ax.plot([dd[L][0] for L in Lo], Lo, "o", mfc="w", mec="k", ms=4)
    ax.set_xticks([1e-3, 3e-3, 1e-2]); ax.set_xticklabels(["$10^{-3}$", r"$3\!\times\!10^{-3}$", "$10^{-2}$"], fontsize=7)
    if k > 0: ax.set_ylabel(""); ax.set_yticklabels([])
fig.text(0.02, 0.96, "(a)", fontsize=10, weight="bold")
ax = fig.add_subplot(gs[1]); Lline2 = np.linspace(9.5, 18.5, 60)
for k, th in enumerate(TH_PROD):
    L_, D_, lab_, ff_, _ = regime_table(fixed, th, T0); dd = {L: crossing(D_, ff_[i]) for i, L in enumerate(L_)}
    Lc = [L for L in L_ if dd[L][1] == 0]
    if len(Lc) >= 2:
        ax.plot(Lc, [dd[L][0] for L in Lc], "o", ms=4, color=f"C{k}", label=rf"$\theta = {th}$")
        Lr = np.linspace(min(Lc), max(Lc), 40); ax.plot(Lr, [free.Dc(L, th=th) for L in Lr], "-", lw=1.1, color=f"C{k}")
ax.plot(Lline2, Lline2 ** 2 / (np.pi ** 2 * T0), "k:", lw=1.3, label=r"$L^2/(\pi^2 T)$")
ax.plot([], [], "-", color="gray", lw=1.1, label="free-space reference")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Domain size $L$"); ax.set_ylabel(r"Crossover boundary $D_c(L;\theta)$")
ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.set_ylim(2.2e-4, 1.2e-1); ax.legend(ncol=4, fontsize=7.5, loc="lower center"); ax.grid(alpha=0.3, which="both")
fig.text(0.02, 0.45, "(b)", fontsize=10, weight="bold")
fig.legend(handles=[Patch(color=RCOL[r], label=r) for r in REGIMES], loc="upper right", ncol=4, fontsize=7.5, frameon=False, bbox_to_anchor=(0.98, 0.995))
save(fig, "threshold_sensitivity")
# agreement of the free-space curves with the boundaries at the other thresholds
relth = []
for th in TH_PROD:
    L_, D_, lab_, ff_, _ = regime_table(fixed, th, T0)
    for i, L in enumerate(L_):
        dcc = crossing(D_, ff_[i])
        if dcc[1] == 0:
            dfree = free.Dc(L, th=th, grid=D_)
            if not np.isnan(dfree): relth.append((th, L, 100 * (dcc[0] / dfree - 1)))
relth = np.array(relth); print("boundary vs free-space, all thresholds: max |rel| per theta:", {th: round(np.max(np.abs(relth[relth[:, 0] == th, 2])), 1) for th in TH_PROD if np.any(relth[:, 0] == th)})
N("bdyRelAllTheta", up(np.max(np.abs(relth[relth[:, 0] <= 0.05, 2])), 0)); N("bdyRelTen", up(np.max(np.abs(relth[relth[:, 0] == 0.1, 2])), 0)); N("bdyRelTwenty", f"{np.max(np.abs(relth[relth[:, 0] == 0.2, 2])):.0f}")
N("LstarTen", f"{np.sqrt(2 * free.A_grid(T0, 0.1).max()):.1f}"); N("LstarTwenty", f"{np.sqrt(2 * free.A_grid(T0, 0.2).max()):.1f}")
print("theta = 0.1, 0.2 boundary points (theta, L, rel%):", [tuple(np.round(r_, 2)) for r_ in relth if r_[0] >= 0.1], "| L* for theta = 0.1, 0.2:", [round(float(np.sqrt(2 * free.A_grid(T0, th_).max())), 1) for th_ in (0.1, 0.2)])

# ---------- Fig. (hotspot size and horizon) ----------
variants = [("r0_0.25", 0.25, 1.0, "v", "C0"), ("eps_0.5", 0.5, 0.5, "d", "C2"), ("baseline", 0.5, 1.0, "o", "k"), ("r0_1.0", 1.0, 1.0, "^", "C3")]
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6)); ax = axs[0]; hsrel = []
for tag, r0, eps, mk, col in variants:
    runs = hi if tag == "baseline" else load_group(os.path.join(R1, "hotsize", tag))
    fr = free if tag == "baseline" else Free(tag)
    L_, D_, lab_, ff_, _ = regime_table(runs, 0.01, T0); dd = {L: crossing(D_, ff_[i]) for i, L in enumerate(L_)}
    Lc = np.array([L for L in L_ if dd[L][1] == 0]); Dc_ = np.array([dd[L][0] for L in Lc])
    c_, l0_, e_ = offsetfit(Lc, Dc_)
    ax.plot(Lc, Dc_, mk, color=col, ms=5, label=rf"$r_0 = {r0:g}$, $\epsilon = {eps:g}$ ($L_0 = {l0_:.1f}$)")
    if len(fr.Dg):
        Lr = np.linspace(Lc.min(), Lc.max(), 50); ax.plot(Lr, [fr.Dc(L) for L in Lr], "-", color=col, lw=1.1)
        hsrel += [100 * abs(d / fr.Dc(L, grid=D_) - 1) for L, d in zip(Lc, Dc_) if not np.isnan(fr.Dc(L, grid=D_))]
ax.plot([], [], "-", color="gray", lw=1.1, label="free-space reference")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.set_xlabel("Domain size $L$"); ax.set_ylabel("$D_c(L)$")
ax.set_title("(a) Hotspot size, $T = 400$", fontsize=9); ax.legend(fontsize=6.3); ax.grid(alpha=0.3, which="both")
if hsrel: N("hsRel", up(max(hsrel), 1)); print(f"hotspot variants: boundary vs free-space max rel {max(hsrel):.1f}%")
ax = axs[1]; local = load_group(os.path.join(R1, "local"))
for T, mk in ((200.0, "s"), (300.0, "D"), (400.0, "o")):
    L_, D_, lab_, ff_, _ = regime_table(fixed, 0.01, T); dd = {L: crossing(D_, ff_[i]) for i, L in enumerate(L_)}
    Lc = np.array([L for L in L_ if dd[L][1] == 0]); ax.plot(Lc, np.array([dd[L][0] for L in Lc]) * T, mk, color="k", ms=4, mfc="none" if T != 400 else "k", label=rf"hotspot, $T' = {T:g}$")
    Lr = np.linspace(Lc.min(), Lc.max(), 40); ax.plot(Lr, np.array([free.Dc(L, T=T) for L in Lr]) * T, "-", color="k", lw=0.8)
for T, mk in ((15.0, "s"), (20.0, "D"), (25.0, "o")):
    L_, D_, lab_, ff_, _ = regime_table(local, 0.01, T); dd = {L: crossing(D_, ff_[i]) for i, L in enumerate(L_)}
    Lc = np.array([L for L in L_ if dd[L][1] == 0]); ax.plot(Lc, np.array([dd[L][0] for L in Lc]) * T * T, mk, color="C0", ms=4, mfc="none" if T != 25 else "C0", label=rf"local, $T' = {T:g}$")
ax.plot(Lline2, Lline2 ** 2 / (8 * np.pi), "--", color="C0", lw=1, label=r"$L^2/(8\pi k_1 S_b)$")
ax.plot([], [], "-", color="k", lw=0.8, label="free-space reference")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.set_xlabel("Domain size $L$"); ax.set_ylabel(r"$D_c\,T'$ (hotspot),  $D_c\,T'^2$ (local)")
ax.set_title("(b) Horizon dependence", fontsize=9); ax.legend(fontsize=6.3, ncol=2); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); save(fig, "hotspot_horizon")

# ---------- Fig. (peak area): add the free-space coverage function ----------
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 7.2))
for i, L in enumerate(Ls):
    a1.plot(Ds, pk[i], "o-", ms=3.5, lw=1, label=f"$L={L:g}$"); a2.plot(Ds, pk[i], "o-", ms=3.5, lw=1, label=f"$L={L:g}$")
p18, c18 = powerfit(Ds, pk[-1]); xx = np.array([Ds[0], Ds[-1]])
Dr = np.geomspace(Ds[0], Ds[-1], 80)
for a in (a1, a2): a.plot(Dr, [free.A(D) for D in Dr], "-", color="k", lw=1.6, label=r"free space, $A_{\rm free}(T;D_I)$")
a2.plot(xx, c18 * xx ** p18, "k--", lw=1.2, label=rf"fit $L={Ls[-1]:g}$: $D_I^{{{p18:.2f}}}$")
a1.set_xlabel("Diffusion coefficient $D_I$"); a1.set_ylabel(r"Peak infected area $A_I^{\rm peak}$"); a1.set_title("(a) Diffusion-dependent peak infected area", fontsize=9)
a2.set_xscale("log"); a2.set_yscale("log"); a2.set_xlabel("Diffusion coefficient $D_I$"); a2.set_ylabel(r"$A_I^{\rm peak}$"); a2.set_title("(b) Log-log scaling", fontsize=9)
a1.legend(ncol=3, fontsize=7); a2.legend(ncol=3, fontsize=7); a1.grid(alpha=0.3); a2.grid(alpha=0.3, which="both")
fig.tight_layout(); save(fig, "peak_area")

# ---------- Fig. (observables): redrawn from the first-revision run with the legend moved off the curves ----------
obs = load_group(os.path.join(R1, "fixed_obs"))[(15.0, 0.01)]; t_ = obs.t()
fig, axs = plt.subplots(3, 1, figsize=(6.2, 7.3), sharex=True); cm_ = plt.get_cmap("viridis")
for k, th in enumerate(obs.thetas):
    if th in (0.005, 0.02, 0.8): continue
    axs[0].plot(t_, obs.area(th), lw=1.4, color=cm_(k / len(obs.thetas)), label=rf"$\theta = {th:g}$")
axs[0].set_ylabel(r"$A_I(t;\theta)$"); axs[0].legend(ncol=3, fontsize=7, loc="upper left"); axs[0].set_title(rf"(a) Thresholded infected area, $(L, D_I) = ({obs.L:g}, {obs.D:g})$", fontsize=9)
axs[1].plot(t_, obs.ts[:, 1], "C2-", lw=1.6); axs[1].set_ylabel("$M_I(t)$", color="C2")
axb = axs[1].twinx(); axb.plot(t_, obs.ts[:, 3], "C3-", lw=1.6); axb.set_ylabel(r"$I_{\max}(t)$", color="C3")
for th in (0.5, 0.6, 0.7): axb.axhline(th, color="gray", ls=":", lw=0.8)
axs[1].set_title(r"(b) Infected mass and maximum density (dotted: $\theta = 0.5, 0.6, 0.7$)", fontsize=9)
axs[2].plot(t_, obs.ts[:, 5], "C4-", lw=1.6); axs[2].set_ylabel("Participation ratio", color="C4")
axc = axs[2].twinx(); axc.plot(t_, obs.ts[:, 6], "C5-", lw=1.6); axc.set_yscale("log"); axc.set_ylabel("Spatial variance", color="C5")
axs[2].set_title("(c) Uniformity: participation ratio and spatial variance", fontsize=9); axs[2].set_xlabel("Time $t$")
fig.tight_layout(); save(fig, "observables")

# ==========================================================================================================
# extended sweep
# ==========================================================================================================
ext = load_group(os.path.join(DATA, "ext"), need_T=T0)
needed = sorted(tuple(p) for p in json.load(open(os.path.join(HERE, "ext_needed.json"))))
print(f"extended-sweep runs loaded: {len(ext)} of {len(needed)}")
EXT_OK = all((round(float(p[0]), 3), round(float(p[1]), 6)) in ext for p in needed)        # every pair listed in ext_needed.json must have been simulated; only the other cells use the free-space reference
if EXT_OK:
    Le = np.array(sorted({float(p[0]) for p in needed})); De = np.array([round(0.001 * 2 ** (k / 2), 6) for k in range(20)])
    shp = (len(Le), len(De)); cov = np.full(shp, np.nan); covf = cov.copy(); mass = cov.copy(); PR = cov.copy(); src = np.zeros(shp, int)
    labx = np.full(shp, "", dtype=object); drift = []
    for i, L in enumerate(Le):
        for j, D in enumerate(De):
            ts_f = free.ts[D]; Af = ts_f[:, 6 + TH.index(0.01)]
            fr_free_t = np.array([disk_in_square(a, L) for a in Af]) / L ** 2
            covf[i, j] = fr_free_t[-1]
            if (L, D) in ext:
                r = ext[(L, D)]; cov[i, j] = r.area(0.01)[-1] / L ** 2; mass[i, j] = r.ts[-1, 1]; PR[i, j] = r.ts[-1, 5]; src[i, j] = 1
                labx[i, j] = classify(*r.fractions(0.01, T0))
                drift.append(np.max(np.abs(r.ts[:, 1] + r.ts[:, 2] - r.ts[0, 1] - r.ts[0, 2])) / (r.ts[0, 1] + r.ts[0, 2]))
            else:                                                # walls not felt: the confined run equals the free-space reference
                cov[i, j] = covf[i, j]; mass[i, j] = ts_f[-1, 1]; labx[i, j] = classify(fr_free_t.max(), fr_free_t[-1])
                rr, II = free.prof[D]; i2 = np.trapz(2 * np.pi * rr * II * II, rr); PR[i, j] = mass[i, j] ** 2 / (L * L * i2)
    N("extN", len(ext)); N("extDrift", fmt_e(max(drift), 0)); print(f"ext: {len(ext)} runs, max mass drift {max(drift):.1e}")
    # cross-check with the production sweep (even L = 10..18), interpolating the production coverage in log D
    chk = [abs(cov[i, j] - np.interp(np.log(D), np.log(Ds), ff[list(Ls).index(L)])) for i, L in enumerate(Le) if L in Ls for j, D in enumerate(De) if Ds[0] <= D <= Ds[-1] and src[i, j]]
    N("extChk", up(max(chk), 3)); print(f"ext vs production coverage: max |diff| {max(chk):.3f} over {len(chk)} points")
    extchk = load_group(os.path.join(DATA, "extchk"), need_T=T0); ec = []
    for (L, D), r in sorted(extchk.items()):
        if (L, D) in ext:
            a_, b_ = ext[(L, D)].area(0.01)[-1] / L ** 2, r.area(0.01)[-1] / L ** 2; ec.append(abs(a_ - b_)); print(f"   extchk L={L:g} D={D:g}: coverage coarse {a_:.4f} fine {b_:.4f}; PR {ext[(L,D)].ts[-1,5]:.4f} vs {r.ts[-1,5]:.4f}")
    N("extFine", up(max(ec), 4) if ec else "--"); N("extFineN", {0: "no", 1: "one", 2: "two", 3: "three", 4: "four"}.get(len(ec), str(len(ec))))
    # coverage boundary of the extended sweep and comparison with the free-space level curve on the same grid
    dce = {L: crossing(De, cov[i]) for i, L in enumerate(Le)}
    relx = {L: 100 * (dce[L][0] / free.Dc(L, grid=De) - 1) for L in Le if dce[L][1] == 0 and not np.isnan(free.Dc(L, grid=De))}
    print("ext boundary:", {int(L): (round(dce[L][0], 4), dce[L][1]) for L in Le}); print("ext boundary vs free-space (same grid) %:", {int(L): round(v, 1) for L, v in relx.items()})
    covmax = {L: (cov[i].max(), covf[i].max(), De[int(np.argmax(cov[i]))]) for i, L in enumerate(Le)}
    print("max coverage over D (confined, free∩square, D at max):", {int(L): (round(v[0], 3), round(v[1], 3), v[2]) for L, v in covmax.items()})
    satL = [L for L in Le if covmax[L][0] > 0.5]; nosat = [L for L in Le if covmax[L][0] <= 0.5]
    N("extLsatMax", f"{max(satL):g}"); N("extLnoSatMin", f"{min(nosat):g}" if nosat else "--")
    follow = [L for L in sorted(relx) if abs(relx[L]) <= 5.0]; Lfol = max(follow) if follow else np.nan
    N("extRelLast", f"{abs(relx[max(relx)]):.0f}"); N("extLlast", f"{max(relx):g}"); N("extLfollow", f"{Lfol:g}"); N("extRelFollow", up(max(abs(relx[L]) for L in follow), 1))
    Lm = max(satL); im_ = list(Le).index(Lm); idx = np.where(cov[im_] > 0.5)[0]
    N("extSatNote", f"(at $L = {Lm:g}$ the coverage exceeds one half for ${De[idx[0]]:.2f} \\le D_I \\le {De[idx[-1]]:.2f}$ and its maximum over $D_I$ is ${covmax[Lm][0]:.2f}$, against ${covmax[Lm][1]:.2f}$ for the free-space solution in the same square)")
    win = {L: (De[np.where(cov[i] > 0.5)[0][0]], De[np.where(cov[i] > 0.5)[0][-1]], np.where(cov[i] > 0.5)[0][-1] == len(De) - 1) for i, L in enumerate(Le) if L in satL}
    print("saturation window in D per L:", {int(L): v for L, v in win.items()}); closed = [L for L, v in win.items() if not v[2]]
    N("extWinClosedMin", f"{min(closed):g}" if closed else "--")
    dl = cov - covf; N("extWallMax", f"{100*np.nanmax(dl):.0f}")
    iw, jw = np.unravel_index(np.nanargmax(dl), dl.shape); print(f"max wall effect {100*dl[iw,jw]:.0f}% at L={Le[iw]:g}, D={De[jw]:g} (cov {cov[iw,jw]:.2f}, free {covf[iw,jw]:.2f})")
    cnt = {r_: int(np.sum(labx == r_)) for r_ in REGIMES}; print("extended-sweep labels:", cnt)
    for r_, key in (("saturation", "extCntSat"), ("invasion", "extCntInv"), ("transient", "extCntTra"), ("below detection", "extCntBel")): N(key, cnt[r_])
    words = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    ib, jb = np.where(labx == "below detection"); lowD = [(Le[i], De[j]) for i, j in zip(ib, jb) if De[j] < 0.1]; print("below detection at small D:", lowD, "| at large D:", [(Le[i], De[j]) for i, j in zip(ib, jb) if De[j] >= 0.1])
    it, jt = np.where(labx == "transient"); print("transient:", [(Le[i], De[j]) for i, j in zip(it, jt)])
    assert lowD == [(56.0, 0.001)] and all(De[j] >= 0.5 for j in jt)
    N("extCntTraW", words[cnt["transient"]]); N("extCntBelW", words[cnt["below detection"]]); N("extBelLmin", f"{min(Le[i] for i, j in zip(ib, jb) if De[j] >= 0.1):g}"); N("extTraLmin", f"{min(Le[it]):g}")

    # ---------- homogenization crossover D_h(L): PR(T) = level ----------
    def dh_table(level):
        out = {}
        for i, L in enumerate(Le):
            c_ = crossing(De, PR[i], level)
            if c_[1] == 0: out[L] = c_[0]
        return out
    from scipy.optimize import least_squares, brentq
    def offsetfit_rel(L, D):                                     # c (L - L_h)^2 / (pi^2 T), least squares on the relative (logarithmic) deviations
        f = lambda p: np.log(p[0] * (L - p[1]) ** 2 / (np.pi ** 2 * T0)) - np.log(D); r = least_squares(f, [0.8, 3.0])
        return r.x[0], r.x[1], 100 * np.max(np.abs(np.exp(f(r.x)) - 1))
    LFIT = 40.0                                                  # D_h is fitted where PR = level is crossed well below the maximum of PR over D_I
    fits = {}; prmax = PR.max(axis=1)
    for level in (0.8, 0.9, 0.95):
        dh = dh_table(level); Lh_ = np.array(sorted(dh)); Dh_ = np.array([dh[L] for L in Lh_]); m_ = Lh_ <= LFIT
        c_, l0_, e_ = offsetfit_rel(Lh_[m_], Dh_[m_]); pe_all = powerfit(Lh_[m_], Dh_[m_])[0]; big = m_ & (Lh_ >= 20)
        pe_big = powerfit(Lh_[big], Dh_[big])[0] if big.sum() >= 3 else np.nan
        fits[level] = dict(L=Lh_, D=Dh_, m=m_, c=c_, L0=l0_, err=e_, pe=pe_all, pe_big=pe_big, ratio=Dh_ / (Lh_ ** 2 / (np.pi ** 2 * T0)),
                           dev={L: 100 * (D / (c_ * (L - l0_) ** 2 / (np.pi ** 2 * T0)) - 1) for L, D in zip(Lh_[~m_], Dh_[~m_])})
        print(f"D_h at PR={level}: defined for L {Lh_.min():g}..{Lh_.max():g}; fit on L<={Lh_[m_].max():g}: {c_:.2f} (L-{l0_:.2f})^2/(pi^2 T), max rel err {e_:.1f}%; power exponent {pe_all:.2f}, L>=20 {pe_big:.2f}; ratio to L^2/(pi^2T) {fits[level]['ratio'][m_].min():.2f}..{fits[level]['ratio'][m_].max():.2f}; larger L vs fit: { {int(k): round(v) for k, v in fits[level]['dev'].items()} }")
        print("    D_h:", dict(zip(Lh_.astype(int), np.round(Dh_, 4))))
    print("max of PR(T) over D_I per L:", {int(L): round(float(v), 3) for L, v in zip(Le, prmax)})
    f9 = fits[0.9]; m9 = f9["m"]
    N("DhC", f"{f9['c']:.2f}"); N("DhL", f"{f9['L0']:.1f}"); N("DhErr", up(f9['err'], 0)); N("DhExp", f"{f9['pe']:.1f}"); N("DhExpLarge", f"{f9['pe_big']:.1f}")
    N("DhLmin", f"{f9['L'][m9].min():g}"); N("DhLmax", f"{f9['L'][m9].max():g}"); N("DhRatioMin", f"{f9['ratio'][m9][0]:.2f}"); N("DhRatioMax", f"{f9['ratio'][m9][-1]:.2f}")
    N("DhLLow", f"{fits[0.8]['L0']:.1f}"); N("DhLHigh", f"{fits[0.95]['L0']:.1f}"); N("DhCLow", f"{fits[0.8]['c']:.2f}"); N("DhCHigh", f"{fits[0.95]['c']:.2f}")
    N("DhErrLow", up(fits[0.8]['err'], 0)); N("DhErrHigh", up(fits[0.95]['err'], 0)); N("DhLmaxHigh", f"{fits[0.95]['L'][fits[0.95]['m']].max():g}")
    N("DhSlopeLo", f"{2 * f9['L'][m9].max() / (f9['L'][m9].max() - f9['L0']):.1f}"); N("DhSlopeHi", f"{2 * f9['L'][m9].min() / (f9['L'][m9].min() - f9['L0']):.1f}")   # local logarithmic slope of the fit
    N("DhLoverLh", f"{f9['L'][m9].max() / f9['L0']:.0f}"); N("DhLspan", f"{f9['L'][m9].max() / f9['L'][m9].min():.0f}")
    Lnext = [L for L in f9["dev"]]; N("DhLnext", f"{Lnext[0]:g}" if Lnext else "--"); N("DhDevNext", f"{f9['dev'][Lnext[0]]:.0f}" if Lnext else "--")
    Lnone = [L for L, v in zip(Le, prmax) if v < 0.9]; N("DhLnone", f"{min(Lnone):g}"); N("DhPRmaxNone", f"{prmax[list(Le).index(min(Lnone))]:.2f}"); N("DhPRmaxLast", f"{prmax[-1]:.2f}")
    nonmono = [L for i, L in enumerate(Le) if PR[i, -1] < PR[i].max() - 0.005]; N("DhLnonmono", f"{min(nonmono):g}"); print("PR(T) decreases at the largest D_I for L >=", min(nonmono))
    # instantaneous point release at the centre of the Neumann square: 1/PR = (1 + 2 sum_j exp(-8 q j^2))^2 with q = D pi^2 T / L^2
    qpt = {lev: brentq(lambda q: (1 + 2 * sum(np.exp(-8 * q * j * j) for j in range(1, 40))) ** 2 - 1 / lev, 1e-3, 5) for lev in (0.8, 0.9, 0.95)}
    from scipy.special import erf as _erf
    def _pr_gauss(q):                                            # same release without walls, restricted to the square: Gaussian of variance 2 D T per direction
        sig = np.sqrt(2 * q) / np.pi; m1 = _erf(0.5 / (sig * np.sqrt(2))); q1 = _erf(0.5 / sig) / (2 * np.sqrt(np.pi) * sig); return (m1 * m1 / q1) ** 2
    qfree = brentq(lambda q: _pr_gauss(q) - 0.9, 1e-3, 50); N("DhPtFree", f"{qfree:.2f}"); N("DhPtFreeRatio", f"{qfree / qpt[0.9]:.1f}"); print(f"point release without walls, restricted to the square: PR = 0.9 at q = {qfree:.3f} (ratio {qfree/qpt[0.9]:.2f})")
    N("DhPtLow", f"{qpt[0.8]:.2f}"); N("DhPt", f"{qpt[0.9]:.2f}"); N("DhPtHigh", f"{qpt[0.95]:.2f}"); print("point-release estimate of D_h pi^2 T / L^2:", {k: round(v, 3) for k, v in qpt.items()})
    # ---------- free-space control for the homogenization crossover: PR of the free-space solution restricted to the square ----------
    def pr_free(L, D, n=600):
        r_, I_ = free.prof[D]; x_ = (np.arange(n) + 0.5) * L / n - L / 2; X_, Y_ = np.meshgrid(x_, x_, indexing="ij"); F_ = np.interp(np.hypot(X_, Y_), r_, I_, right=0.0)
        return F_.sum() ** 2 / (n * n * (F_ * F_).sum()), F_.sum() * (L / n) ** 2
    PRf = np.zeros(shp); Mfsq = np.zeros(shp)
    for i, L in enumerate(Le):
        for j, D in enumerate(De): PRf[i, j], Mfsq[i, j] = pr_free(L, D)
    assert np.max(np.abs(PR - PRf)[src == 0]) < 2e-3            # consistency of the surrogate cells
    dPR = PR - PRf; iw_, jw_ = np.unravel_index(np.argmax(dPR), dPR.shape); N("prWallMax", f"{dPR.max():.2f}"); print(f"max PR(confined) - PR(free in square): {dPR.max():.3f} at L={Le[iw_]:g}, D={De[jw_]:g}")
    N("massRetMin", f"{100 * (Mfsq / mass).min():.1f}"); N("extNsurr", int(np.sum(src == 0))); N("extNcells", int(src.size)); print(f"mass of the free-space solution inside the square / confined mass: min {(Mfsq / mass).min():.3f}")
    dhf = {}
    for level in (0.8, 0.9, 0.95):
        out_ = {}
        for i, L in enumerate(Le):
            c_ = crossing(De, PRf[i], level)
            if c_[1] == 0: out_[L] = c_[0]
        dhf[level] = out_; fl = fits[level]; common = [L for L in fl["L"] if L in out_]
        rr_ = [out_[L] / fl["D"][list(fl["L"]).index(L)] for L in common]
        print(f"free-space D_h at PR={level}: defined for L <= {max(out_):g}; ratio to the confined D_h: " + ", ".join(f"{L:g}: {v:.2f}" for L, v in zip(common, rr_)))
        fl["free"] = out_; fl["free_ratio"] = dict(zip(common, rr_))
    fr9 = f9["free_ratio"]; small = [v for L, v in fr9.items() if L <= 24]
    N("DhFreeRatioMin", f"{min(small):.1f}"); N("DhFreeRatioMax", f"{max(small):.1f}"); N("DhFreeLmax", f"{max(fr9):g}"); N("DhFreeRatioLast", f"{fr9[max(fr9)]:.1f}")
    N("DhFreeLnone", f"{min(L for L in Le if L not in f9['free']):g}"); N("DhFreeRhoMin", f"{min(f9['free'][L] * np.pi ** 2 * T0 / L ** 2 for L in fr9 if L <= 24):.2f}"); N("DhFreeRhoMax", f"{max(f9['free'][L] * np.pi ** 2 * T0 / L ** 2 for L in fr9 if L <= 24):.2f}")
    lo_ = list(fits[0.8]["free_ratio"].values()); hi_ = list(fits[0.95]["free_ratio"].values())
    N("DhFreeRatioLowMin", f"{min(lo_[:7]):.1f}"); N("DhFreeRatioLowMax", f"{max(lo_[:7]):.1f}"); N("DhFreeRatioHighMin", f"{min(hi_):.1f}"); N("DhFreeRatioHighMax", f"{max(hi_):.1f}"); N("DhFreeLmaxHigh", f"{max(fits[0.95]['free_ratio']):g}")
    N("prFreeMaxNone", f"{PRf[list(Le).index(min(L for L in Le if L not in f9['free']))].max():.2f}")
    # interpolation error of the boundaries: coarse grid versus the fine free-space grid
    ie_c = max(abs(free.Dc(L, grid=Ds) / free.Dc(L) - 1) for L in Lf); N("DcInterp", up(100*ie_c, 0)); print(f"interpolation error of D_c on the production grid (free-space curve, L=11..18): max {100*ie_c:.1f}%")
    ie_h = []
    for L in (10.0, 14.0, 20.0, 24.0):
        Dfine = [D for D in free.Dg if De[0] <= D <= De[-1]]; prf_f = np.array([pr_free(L, D, n=300)[0] for D in Dfine]); prf_c = np.array([pr_free(L, D, n=300)[0] for D in De])
        ie_h.append(abs(crossing(De, prf_c, 0.9)[0] / crossing(Dfine, prf_f, 0.9)[0] - 1))
    N("DhInterp", up(100*max(ie_h), 0)); print("interpolation error of D_h on the sqrt(2) grid (free-space PR, L=10,14,20,24):", np.round(100 * np.array(ie_h), 2))
    json.dump({"Le": Le.tolist(), "De": De.tolist(), "PR_free_in_square": PRf.tolist(), "mass_free_in_square": Mfsq.tolist()}, open(os.path.join(SUP, "extended_sweep_free_PR.json"), "w"), indent=0)
    rat = [f9['D'][list(f9['L']).index(L)] / dc_conf[L][0] for L in (12.0, 14.0, 16.0, 18.0) if L in f9['L'] and dc_conf[L][1] == 0]
    N("DhOverDcMin", f"{min(rat):.0f}"); N("DhOverDcMax", f"{max(rat):.0f}"); print("D_h / D_c at L = 12..18:", np.round(rat, 1))
    # all-or-nothing coverage in the homogenized range
    hom = PR > 0.95; print("homogenized runs (PR>0.95):", int(hom.sum()), "| coverage values there within 0.02 of 0 or 1:", int(np.sum(hom & ((cov < 0.02) | (cov > 0.98)))))
    N("homN", int(hom.sum())); N("homAllOrNothing", int(np.sum(hom & ((cov < 0.02) | (cov > 0.98)))))
    json.dump({"Le": Le.tolist(), "De": De.tolist(), "coverage": cov.tolist(), "coverage_free_in_square": covf.tolist(), "mass": mass.tolist(), "PR": PR.tolist(),
               "simulated": src.tolist(), "label": labx.tolist()}, open(os.path.join(SUP, "extended_sweep_tables.json"), "w"), indent=0)

    # ---------- new figure: free-space reference and wall effects ----------
    fig, axs = plt.subplots(2, 2, figsize=(7.5, 7.0))
    ax = axs[0, 0]; Dr = np.geomspace(free.Dg[0], free.Dg[-1], 300)
    ax.plot(Dr, [free.A(D) for D in Dr], "k-", lw=1.8, label=r"free space $A_{\rm free}(T;D_I)$")
    for L, col, mk in ((10.0, "C0", "o"), (14.0, "C1", "s"), (18.0, "C2", "^")):
        ip = list(Ls).index(L); ax.plot(Ds, ff[ip] * L * L, mk, color=col, ms=4, mfc="none", label=rf"confined, $L = {L:g}$"); ax.axhline(L * L / 2, color=col, ls=":", lw=0.9)
    for L, col, mk in ((28.0, "C4", "v"), (40.0, "C3", "D")):
        ie = list(Le).index(L); m = src[ie] == 1; ax.plot(De[m], cov[ie][m] * L * L, mk, color=col, ms=4, mfc="none", label=rf"confined, $L = {L:g}$"); ax.axhline(L * L / 2, color=col, ls=":", lw=0.9)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(30, 3000); ax.set_xlabel("Diffusion coefficient $D_I$"); ax.set_ylabel(r"Area above threshold $A_I(T;0.01)$")
    ax.set_title("(a) Coverage function", fontsize=9); ax.legend(fontsize=6.3, loc="upper left"); ax.grid(alpha=0.3, which="both")
    ax = axs[0, 1]
    Lr = np.linspace(10.2, Lstar - 0.01, 400); ax.plot(Lr, [free.Dc(L) for L in Lr], "-", color="C3", lw=1.8, label="free-space level curve")
    ax.plot(Lf, Df, "ko", ms=4, label="confined, production sweep")
    Lin = [L for L in Le if dce[L][1] == 0]; ax.plot(Lin, [dce[L][0] for L in Lin], "s", color="C0", ms=5.5, mfc="none", label="confined, extended sweep")
    Lq = np.linspace(10, 58, 100); ax.plot(Lq, c_off * (Lq - L0) ** 2, "k-", lw=0.9, label=rf"$c\,(L-L_0)^2$, fitted on $L\leq 18$")
    ax.plot(Lq, Lq ** 2 / (np.pi ** 2 * T0), "k:", lw=1.2, label=r"$L^2/(\pi^2 T)$")
    ax.axvline(Lstar, color="C3", ls="--", lw=0.9); ax.text(Lstar + 0.6, 0.6, r"$L^\ast$", color="C3", fontsize=9)
    ax.set_yscale("log"); ax.set_xlim(9, 58); ax.set_ylim(6e-4, 1.5); ax.set_xlabel("Domain size $L$"); ax.set_ylabel("Coverage boundary $D_c(L)$")
    ax.set_title("(b) Coverage boundary, extended range", fontsize=9); ax.legend(fontsize=6.3, loc="lower right"); ax.grid(alpha=0.3, which="both")
    ax = axs[1, 0]; draw_map(ax, Le, De, labx); ax.set_yticks(Le[::2])
    ax.plot([free.Dc(L) for L in Lr], Lr, "-", color="C3", lw=1.6, label="free-space level curve")
    ax.plot([dce[L][0] for L in Lin], Lin, "ko", ms=3.5, label="$D_c$ (confined)")
    ax.plot(f9["D"], f9["L"], "^", color="k", mfc="w", ms=5, label=r"$D_h$: PR$(T) = 0.9$")
    ax.plot([f9["free"][L] for L in sorted(f9["free"])], sorted(f9["free"]), "v", color="C3", mfc="w", ms=5, label=r"PR$(T) = 0.9$, free space")
    Lh = np.linspace(9.5, LFIT, 100); ax.plot(f9["c"] * (Lh - f9["L0"]) ** 2 / (np.pi ** 2 * T0), Lh, "k--", lw=1, label=r"fit $c\,(L-L_h)^2/(\pi^2T)$, $L \leq 40$")
    ii, jj = np.where(dl > 0.01); ax.plot(De[jj], Le[ii], "+", color="k", ms=4, mew=0.7, label="wall effect $> 1\\%$")
    ax.set_xlim(De[0] / 1.2, De[-1] * 1.2); ax.set_ylim(Le[0] - 1, Le[-1] + 2); ax.set_title("(c) Extended regime map", fontsize=9)
    h, l = ax.get_legend_handles_labels(); ax.legend(handles=h + [Patch(color=RCOL[r_], label=r_) for r_ in REGIMES if cnt[r_] > 0], fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
    ax = axs[1, 1]
    for L, col in ((10.0, "C0"), (20.0, "C1"), (32.0, "C2"), (36.0, "C3"), (40.0, "C4"), (48.0, "C5")):
        if L not in Le: continue
        ie = list(Le).index(L); ax.plot(De, cov[ie], "o-", color=col, ms=3, lw=1, label=rf"$L = {L:g}$"); ax.plot(De, covf[ie], "--", color=col, lw=0.9)
    ax.axhline(0.5, color="k", ls=":", lw=0.9); ax.plot([], [], "k--", lw=0.9, label="free space")
    ax.set_xscale("log"); ax.set_xlabel("Diffusion coefficient $D_I$"); ax.set_ylabel(r"Coverage $a_T = A_I(T;0.01)/L^2$"); ax.set_title("(d) Coverage: confined versus free space", fontsize=9)
    ax.legend(fontsize=6.3, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); save(fig, "free_space_walls")

# ----------------------------------------------------------------------------------------------------------
with open(os.path.join(MS, "numbers_r2.tex"), "w") as f:
    f.write("% numbers of the second revision, generated by code/analyze_r2.py -- do not edit by hand\n")
    for k, v in numbers.items(): f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print("numbers:", numbers)
