#!/usr/bin/env python3
"""Classify regimes, fit scalings and produce all figures / numbers for the revision.

Reads  ../data/{fixed,fixed_hi_theta,fixed_obs,local,kernel,conv,long}/.../{summary.txt,timeseries.txt,snap_*.bin}
Writes ../manuscript/fig2.eps ... fig8.eps (+png), numbers.tex, table_conv.tex, ../data/regimes_*.csv
"""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
MS = os.path.normpath(os.path.join(HERE, "..", "manuscript"))
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 150, "savefig.dpi": 300})

# ----------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------
def read_timeseries(path):
    thetas, meta = None, {}
    with open(path) as f:
        for line in f:
            if line.startswith("# thetas:"): thetas = np.array([float(x) for x in line.split(":")[1].split()])
            elif line.startswith("# model="):
                for kv in line[2:].split(): k, v = kv.split("="); meta[k] = v
            elif not line.startswith("#"): break
    arr = np.loadtxt(path, comments="#")
    cols = {"t": 0, "massI": 1, "massS": 2, "maxI": 3, "centerI": 4, "PR": 5, "var": 6}
    return arr, cols, thetas, meta

class Run:
    def __init__(self, d):
        self.dir = d
        self.ts, self.cols, self.thetas, self.meta = read_timeseries(os.path.join(d, "timeseries.txt"))
        self.L = float(self.meta["L"]); self.D = float(self.meta["D"]); self.model = self.meta["model"]
        self.Nx = int(self.meta["Nx"]); self.dx = float(self.meta["dx"])
    def area(self, theta): return self.ts[:, 7 + int(np.argmin(np.abs(self.thetas - theta)))]
    def t(self): return self.ts[:, 0]
    def fractions(self, theta, T):
        """(peak fraction over [0,T], fraction at T)"""
        a = self.area(theta) / self.L ** 2; j = int(np.argmin(np.abs(self.t() - T)))
        return a[: j + 1].max(), a[j]

def load_group(group):
    runs = {}
    for d in sorted(glob.glob(os.path.join(DATA, group, "*L_*_D_*"))):
        if not os.path.exists(os.path.join(d, "summary.txt")): continue
        r = Run(d); runs[(round(r.L, 3), round(r.D, 6))] = r
    return runs

# ----------------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------------
def classify(peak_frac, final_frac):
    if peak_frac < 0.02: return "collapse"
    if final_frac > 0.50: return "saturation"
    if final_frac < 0.5 * peak_frac: return "transient"
    return "invasion"

REGIMES = ["collapse", "transient", "invasion", "saturation"]
RCOL = {"collapse": "#d9d9d9", "transient": "#fdb462", "invasion": "#80b1d3", "saturation": "#8dd3c7"}

def regime_table(runs, theta, T):
    Ls = sorted({k[0] for k in runs}); Ds = sorted({k[1] for k in runs})
    lab = np.full((len(Ls), len(Ds)), "", dtype=object)
    ff = np.full((len(Ls), len(Ds)), np.nan); pk = np.full((len(Ls), len(Ds)), np.nan)
    for (L, D), r in runs.items():
        pf, fr = r.fractions(theta, T); i, j = Ls.index(L), Ds.index(D)
        lab[i, j] = classify(pf, fr); ff[i, j] = fr; pk[i, j] = pf * L * L
    return np.array(Ls), np.array(Ds), lab, ff, pk

def crossover(Ls, Ds, ff, level=0.5):
    """interpolated D_c(L) where the final fraction crosses `level` (linear in log D).
    censor flag: -1 below grid, +1 above grid, 0 inside."""
    out = {}
    for i, L in enumerate(Ls):
        f = ff[i]; ok = ~np.isnan(f); d, f = Ds[ok], f[ok]
        if len(d) == 0: continue
        if f[0] > level: out[L] = (d[0], -1); continue
        if f[-1] <= level: out[L] = (d[-1], +1); continue
        j = int(np.argmax(f > level))
        x0, x1, y0, y1 = np.log(d[j - 1]), np.log(d[j]), f[j - 1], f[j]
        out[L] = (float(np.exp(x0 + (level - y0) * (x1 - x0) / (y1 - y0))), 0)
    return out

def inside(dc):
    Lf = np.array([L for L in dc if dc[L][1] == 0]); return Lf, np.array([dc[L][0] for L in Lf])

def powerfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); m = (x > 0) & (y > 0)
    lx, ly = np.log(x[m]), np.log(y[m]); p, b = np.polyfit(lx, ly, 1); pred = p * lx + b
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2) if len(lx) > 2 else np.nan
    return p, np.exp(b), r2

def offsetfit(L, D):
    f = lambda L, c, L0: c * (L - L0) ** 2
    popt, _ = curve_fit(f, L, D, p0=[1e-5, 5.0]); pred = f(L, *popt)
    return popt[0], popt[1], 100 * np.max(np.abs(pred - D) / D)

def draw_map(ax, Ls, Ds, lab, title=None, dots=None):
    Z = np.array([[REGIMES.index(l) if l else np.nan for l in row] for row in lab], float)
    lx = np.log10(Ds); e = np.zeros(len(Ds) + 1)
    e[1:-1] = 0.5 * (lx[:-1] + lx[1:]); e[0] = lx[0] - (e[1] - lx[0]); e[-1] = lx[-1] + (lx[-1] - e[-2])
    ye = np.zeros(len(Ls) + 1); ye[1:-1] = 0.5 * (Ls[:-1] + Ls[1:]); ye[0] = Ls[0] - (ye[1] - Ls[0]); ye[-1] = Ls[-1] + (Ls[-1] - ye[-2])
    cmap = ListedColormap([RCOL[r] for r in REGIMES]); norm = BoundaryNorm(np.arange(-0.5, 4.5, 1), 4)
    ax.pcolormesh(10 ** e, ye, Z, cmap=cmap, norm=norm, shading="flat", edgecolors="w", linewidth=0.3)
    ax.set_xscale("log"); ax.set_xlabel("Diffusion coefficient $D_I$"); ax.set_ylabel("Domain size $L$"); ax.set_yticks(Ls)
    if title: ax.set_title(title, fontsize=9)
    if dots:
        Lc, Dc = inside(dots)
        if len(Lc): ax.plot(Dc, Lc, "ko", ms=4, label="numerical boundary $D_c(L)$")
        Lc2 = [L for L in dots if dots[L][1] != 0]
        if Lc2: ax.plot([dots[L][0] for L in Lc2], Lc2, "o", mfc="w", mec="k", ms=4, label="boundary outside grid")

def fmt_e(x):
    m, e = f"{x:.1e}".split("e"); return rf"{m}\times 10^{{{int(e)}}}"

numbers = {}
def N(key, val): numbers[key] = val

# ============================================================================
fixed = load_group("fixed"); T_FIXED = 400.0; print(f"fixed runs: {len(fixed)}")
Ls, Ds, lab, ff, pk = regime_table(fixed, 0.01, T_FIXED)
cnt = {r: int(np.sum(lab == r)) for r in REGIMES}; print("theta=0.01 counts:", cnt)
N("nInv", cnt["invasion"]); N("nSat", cnt["saturation"])
dc = crossover(Ls, Ds, ff); Lf, Df = inside(dc)
pexp, pref, r2 = powerfit(Lf, Df); c_off, L0, err_off = offsetfit(Lf, Df)
pure = 1.0 / (np.pi ** 2 * T_FIXED)
print(f"D_c: power exp {pexp:.2f} (R2 {r2:.3f}); offset fit c={c_off:.3e} L0={L0:.2f} maxerr {err_off:.1f}%; 1/(pi^2 T)={pure:.2e}")
N("DcExp", f"{pexp:.1f}"); N("DcOffC", fmt_e(c_off)); N("DcOffL", f"{L0:.1f}"); N("DcOffErr", f"{err_off:.0f}")
N("DcOffRatio", f"{pure / c_off:.0f}"); N("DcMin", f"{Df.min():.4f}"); N("DcMax", f"{Df.max():.4f}")
N("DcLmin", f"{Lf.min():g}"); N("DcLmax", f"{Lf.max():g}")
# hotspot radius where chi = 0.05 and mass produced
r05 = 0.5 + np.arctanh(1 - 2 * 0.05); N("hotR", f"{r05:.1f}")
mI = np.array([[fixed[(L, D)].ts[-1, 1] for D in Ds] for L in Ls]); N("massMin", f"{mI.min():.0f}"); N("massMax", f"{mI.max():.0f}")

# ---------------- Fig 2 ----------------
fig, ax = plt.subplots(figsize=(6.4, 4.4))
draw_map(ax, Ls, Ds, lab, dots=dc)
Lline = np.linspace(Ls[0] - 0.5, Ls[-1] + 0.5, 60)
ax.plot(c_off * (Lline - L0) ** 2, Lline, "k-", lw=1.4, label=rf"$D_c = c\,(L - L_0)^2$, $L_0 = {L0:.1f}$")
ax.plot(pref * Lline ** pexp, Lline, "k--", lw=1.2, label=rf"$D_c \propto L^{{{pexp:.1f}}}$")
ax.set_xlim(10 ** (np.log10(Ds[0]) - 0.06), 10 ** (np.log10(Ds[-1]) + 0.06)); ax.set_ylim(Ls[0] - 0.5, Ls[-1] + 0.5)
h, l = ax.get_legend_handles_labels(); h += [Patch(color=RCOL["invasion"], label="invasion"), Patch(color=RCOL["saturation"], label="saturation")]
ax.legend(handles=h, loc="upper left", fontsize=7.5); ax.set_title(r"Regime map, hotspot model ($\theta = 0.01$, $T = 400$)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig2.eps")); fig.savefig(os.path.join(MS, "fig2.png")); plt.close(fig)

# ---------------- Fig 3 ----------------
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 7.2)); alphas = {}
for i, L in enumerate(Ls):
    y = pk[i]; m = ~np.isnan(y)
    a1.plot(Ds[m], y[m], "o-", ms=3.5, lw=1, label=f"$L={L:g}$"); a2.plot(Ds[m], y[m], "o-", ms=3.5, lw=1, label=f"$L={L:g}$")
    if L >= 16: alphas[L] = powerfit(Ds[m], y[m])
aL = Ls[-1]; p, c, r = alphas[aL]
N("alphafit", f"{p:.2f}"); N("alphaR", f"{r:.3f}"); N("alphaSixteen", f"{alphas[16.0][0]:.2f}"); N("alphaSeventeen", f"{alphas[17.0][0]:.2f}")
print("alpha:", {k: round(v[0], 3) for k, v in alphas.items()})
xx = np.array([Ds[0], Ds[-1]]); a2.plot(xx, c * xx ** p, "k--", lw=1.5, label=rf"fit $L={aL:g}$: $D_I^{{{p:.2f}}}$")
a1.set_xlabel("Diffusion coefficient $D_I$"); a1.set_ylabel(r"Peak infected area $A_I^{\rm peak}$"); a1.set_title("(a) Diffusion-dependent peak infected area", fontsize=9)
a2.set_xscale("log"); a2.set_yscale("log"); a2.set_xlabel("Diffusion coefficient $D_I$"); a2.set_ylabel(r"$A_I^{\rm peak}$"); a2.set_title("(b) Log-log scaling", fontsize=9)
a1.legend(ncol=3, fontsize=7); a2.legend(ncol=3, fontsize=7); a1.grid(alpha=0.3); a2.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig5.eps")); fig.savefig(os.path.join(MS, "fig5.png")); plt.close(fig)

# ---------------- Fig 4 ----------------
def load_snap(r, t): return np.fromfile(os.path.join(r.dir, f"snap_t{t:08.2f}.bin"), dtype=np.float64).reshape(r.Nx, r.Nx)
rep = fixed[(15.0, 0.01)]
fig = plt.figure(figsize=(7.2, 6.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], left=0.08, right=0.87, top=0.95, bottom=0.08, hspace=0.5, wspace=0.35)
snaps = [load_snap(rep, t) for t in (100, 200, 400)]; vmax = max(s.max() for s in snaps)
for k, (t, sn) in enumerate(zip((100, 200, 400), snaps)):
    ax = fig.add_subplot(gs[0, k]); im = ax.imshow(sn.T, origin="lower", extent=[0, rep.L, 0, rep.L], vmin=0, vmax=vmax, cmap="viridis")
    ax.set_title(f"({'abc'[k]}) $t = {t}$", fontsize=9); ax.set_xlabel("$x$")
    if k == 0: ax.set_ylabel("$y$")
cax = fig.add_axes([0.89, 0.60, 0.015, 0.3]); cb = fig.colorbar(im, cax=cax); cb.set_label("Infected density $I$")
ax = fig.add_subplot(gs[1, :]); t = rep.t()
ax.plot(t, rep.area(0.01), "C0-", lw=1.8, label=r"$A_I(t;0.01)$"); ax.set_xlabel("Time $t$"); ax.set_ylabel("Infected area $A_I$", color="C0"); ax.set_title("(d) Macroscopic observables", fontsize=9)
ax2 = ax.twinx(); ax2.plot(t, rep.ts[:, 1], "C2--", lw=1.5, label="$M_I(t)$"); ax2.plot(t, rep.ts[:, 2] / 10, "C1-.", lw=1.5, label="$M_S(t)/10$")
ax2.plot(t, rep.ts[:, 3] * 100, "C3:", lw=1.8, label=r"$100\,I_{\max}(t)$"); ax2.set_ylabel("Mass / scaled maximum density")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, loc="center right", fontsize=7.5)
fig.savefig(os.path.join(MS, "fig6.eps")); fig.savefig(os.path.join(MS, "fig6.png")); plt.close(fig)
N("repImaxPeak", f"{rep.ts[:, 3].max():.2f}"); N("repImaxPeakT", f"{t[np.argmax(rep.ts[:, 3])]:.0f}"); N("repImaxT", f"{rep.ts[-1, 3]:.2f}")

# ---------------- Fig 5 ----------------
H = np.array([fixed[(L, 0.01)].area(0.01) for L in Ls]); tgrid = fixed[(Ls[0], 0.01)].t()
fig, ax = plt.subplots(figsize=(6.4, 4.4)); im = ax.pcolormesh(tgrid, Ls, H, shading="nearest", cmap="plasma")
t50 = [tgrid[np.argmax(H[i] / L ** 2 >= 0.5)] if np.any(H[i] / L ** 2 >= 0.5) else np.nan for i, L in enumerate(Ls)]
ax.plot(t50, Ls, "w-", lw=2, label="area fraction $= 0.5$"); ax.plot(t50, Ls, "wo", ms=3)
cb = fig.colorbar(im, ax=ax); cb.set_label("Infected area $A_I(t;0.01)$"); ax.set_xlabel("Time $t$"); ax.set_ylabel("Domain size $L$"); ax.set_yticks(Ls); ax.legend(loc="lower right", fontsize=8)
ax.set_title(r"Infected area across domain sizes, $D_I = 0.01$", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig8.eps")); fig.savefig(os.path.join(MS, "fig8.png")); plt.close(fig)

# ---------------- Fig 6: threshold sensitivity ----------------
hi = load_group("fixed_hi_theta"); print(f"coarse runs: {len(hi)}")
TH_PROD = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]; TH_HI = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
dcs, exps, cnts = {}, {}, {}
for th in TH_PROD:
    L_, D_, lab_, ff_, _ = regime_table(fixed, th, T_FIXED); dcs[th] = crossover(L_, D_, ff_)
    Lf_, Df_ = inside(dcs[th]); cnts[th] = {r: int(np.sum(lab_ == r)) for r in REGIMES}
    if len(Lf_) >= 3: exps[th] = powerfit(Lf_, Df_)[0]
    print(f"theta={th}: {cnts[th]}, exponent {exps.get(th, np.nan):.2f}")
cnts_hi = {}
if hi:
    for th in TH_HI:
        L_, D_, lab_, ff_, _ = regime_table(hi, th, T_FIXED); cnts_hi[th] = {r: int(np.sum(lab_ == r)) for r in REGIMES}
        print(f"coarse theta={th}: {cnts_hi[th]}")
    # coarse-vs-production check at theta = 0.01
    L_, D_, labh, ffh, _ = regime_table(hi, 0.01, T_FIXED); same = np.mean(labh == lab); dmax = np.nanmax(np.abs(ffh - ff))
    print(f"coarse vs production: identical labels {100*same:.0f}%, max |delta a_T| = {dmax:.3f}"); N("coarseAgree", f"{100*same:.0f}"); N("coarseDa", f"{dmax:.3f}")
fig = plt.figure(figsize=(7.2, 7.0)); gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.95], hspace=0.45, left=0.1, right=0.97, top=0.9, bottom=0.07)
gsa = gs[0].subgridspec(1, 4, wspace=0.12)
panels = [(0.005, fixed), (0.05, fixed), (0.3, hi if hi else fixed), (0.6, hi)] if hi else [(0.005, fixed), (0.05, fixed), (0.1, fixed), (0.3, fixed)]
for k, (th, grp) in enumerate(panels):
    ax = fig.add_subplot(gsa[0, k]); L_, D_, lab_, ff_, _ = regime_table(grp, th, T_FIXED)
    draw_map(ax, L_, D_, lab_, title=rf"$\theta = {th}$", dots=crossover(L_, D_, ff_))
    ax.set_xticks([1e-3, 3e-3, 1e-2]); ax.set_xticklabels(["$10^{-3}$", r"$3\!\times\!10^{-3}$", "$10^{-2}$"], fontsize=7)
    if k > 0: ax.set_ylabel(""); ax.set_yticklabels([])
fig.text(0.02, 0.96, "(a)", fontsize=10, weight="bold")
ax = fig.add_subplot(gs[1])
for th in TH_PROD:
    Lf_, Df_ = inside(dcs[th])
    if len(Lf_) >= 2: ax.plot(Lf_, Df_, "o-", ms=4, lw=1.2, label=rf"$\theta = {th}$")
ax.plot(Lline, pure * Lline ** 2, "k:", lw=1.3, label=r"$L^2/(\pi^2 T)$")
ax.plot(Lline[Lline > L0 + 1], c_off * (Lline[Lline > L0 + 1] - L0) ** 2, "k-", lw=1, alpha=0.6, label=rf"$c\,(L-L_0)^2$, $\theta = 0.01$")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Domain size $L$"); ax.set_ylabel(r"Crossover boundary $D_c(L;\theta)$")
ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.legend(ncol=3, fontsize=7.5); ax.grid(alpha=0.3, which="both")
fig.text(0.02, 0.45, "(b)", fontsize=10, weight="bold")
fig.legend(handles=[Patch(color=RCOL[r], label=r) for r in REGIMES], loc="upper right", ncol=4, fontsize=7.5, frameon=False, bbox_to_anchor=(0.98, 0.995))
fig.savefig(os.path.join(MS, "fig3.eps")); fig.savefig(os.path.join(MS, "fig3.png")); plt.close(fig)
ratio = []
for L in Ls:
    if all(dcs[t_][L][1] == 0 for t_ in (0.005, 0.02)): ratio.append(dcs[0.02][L][0] / dcs[0.005][L][0])
N("thetaRatio", f"{max(ratio):.2f}" if ratio else "--"); N("DcExpMin", f"{min(exps[t_] for t_ in exps if t_ <= 0.1):.1f}"); N("DcExpMax", f"{max(exps[t_] for t_ in exps if t_ <= 0.1):.1f}")
N("satPointOne", f"{cnts[0.1]['saturation']}"); N("satPointThree", f"{cnts[0.3]['saturation']}")
if hi:
    N("transSix", f"{cnts_hi[0.6]['transient']}"); N("transSeven", f"{cnts_hi[0.7]['transient']}"); N("collSeven", f"{cnts_hi[0.7]['collapse']}"); N("collEight", f"{cnts_hi[0.8]['collapse']}")
    N("transFive", f"{cnts_hi[0.5]['transient']}"); N("satSix", f"{cnts_hi[0.6]['saturation']}"); N("invSix", f"{cnts_hi[0.6]['invasion']}")

# ---------------- Fig 7: observables on one run ----------------
obsg = load_group("fixed_obs"); obs = obsg[(15.0, 0.01)] if (15.0, 0.01) in obsg else rep
t = obs.t(); th_all = obs.thetas
fig, axs = plt.subplots(3, 1, figsize=(6.2, 7.3), sharex=True)
cm = plt.get_cmap("viridis")
for k, th in enumerate(th_all):
    if th in (0.005, 0.02, 0.8): continue
    axs[0].plot(t, obs.area(th), lw=1.4, color=cm(k / len(th_all)), label=rf"$\theta = {th:g}$")
axs[0].set_ylabel(r"$A_I(t;\theta)$"); axs[0].legend(ncol=5, fontsize=7); axs[0].set_title(rf"(a) Thresholded infected area, $(L, D_I) = ({obs.L:g}, {obs.D:g})$", fontsize=9)
axs[1].plot(t, obs.ts[:, 1], "C2-", lw=1.6, label="$M_I(t)$"); axs[1].set_ylabel("$M_I(t)$", color="C2")
axb = axs[1].twinx(); axb.plot(t, obs.ts[:, 3], "C3-", lw=1.6); axb.set_ylabel(r"$I_{\max}(t)$", color="C3")
for th in (0.5, 0.6, 0.7): axb.axhline(th, color="gray", ls=":", lw=0.8)
axs[1].set_title(r"(b) Infected mass and maximum density (dotted: $\theta = 0.5, 0.6, 0.7$)", fontsize=9)
axs[2].plot(t, obs.ts[:, 5], "C4-", lw=1.6); axs[2].set_ylabel("Participation ratio", color="C4")
axc = axs[2].twinx(); axc.plot(t, obs.ts[:, 6], "C5-", lw=1.6); axc.set_yscale("log"); axc.set_ylabel("Spatial variance", color="C5")
axs[2].set_title("(c) Homogenization: participation ratio and spatial variance", fontsize=9); axs[2].set_xlabel("Time $t$")
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig7.eps")); fig.savefig(os.path.join(MS, "fig7.png")); plt.close(fig)
N("obsL", f"{obs.L:g}"); N("obsD", f"{obs.D:g}"); N("obsImax", f"{obs.ts[-1, 3]:.2f}"); N("obsImaxPeak", f"{obs.ts[:, 3].max():.2f}"); N("obsImaxPeakT", f"{t[np.argmax(obs.ts[:, 3])]:.0f}")
N("obsMean", f"{obs.ts[-1, 1] / obs.L ** 2:.3f}"); N("obsPR", f"{obs.ts[-1, 5]:.2f}"); N("obsMass", f"{obs.ts[-1, 1]:.1f}")
N("obsMassHalfT", f"{t[np.argmax(obs.ts[:, 1] >= 0.95 * obs.ts[-1, 1])]:.0f}")
N("obsMassAtPeak", f"{obs.ts[np.argmin(np.abs(t - 150)), 1]:.0f}")
N("obsFracFive", f"{obs.area(0.05)[-1] / obs.L ** 2:.2f}"); N("hotDiam", f"{2 * r05:.1f}"); N("obsMassTwoHundred", f"{obs.ts[np.argmin(np.abs(t - 200)), 1]:.0f}"); N("obsMassFinal", f"{obs.ts[-1, 1]:.0f}")
decay = []
for th in th_all:
    a = obs.area(th); pk_ = a.max()
    if pk_ > 0 and a[-1] < 0.05 * pk_: decay.append((th, pk_, t[np.argmax(a)], a[-1]))
print("rise-and-decay thresholds:", decay)
if decay:
    N("obsDecayPeakFrac", f"{100 * decay[0][1] / obs.L ** 2:.0f}\\%"); N("obsDecayThetas", ", ".join(f"{d[0]:g}" for d in decay)); N("obsDecayFirst", f"{decay[0][0]:g}"); N("obsDecayPeak", f"{decay[0][1]:.1f}"); N("obsDecayT", f"{decay[0][2]:.0f}")
    N("obsDecayEnd", f"{t[np.argmax(obs.area(decay[0][0])[np.argmax(obs.area(decay[0][0])):] < 0.05 * decay[0][1]) + np.argmax(obs.area(decay[0][0]))]:.0f}")
else:
    for k in ("obsDecayThetas", "obsDecayFirst", "obsDecayPeak", "obsDecayT", "obsDecayEnd"): N(k, "--")
# which theta between late and early maximum
N("obsPlateau", ", ".join(f"{th:g}" for th in th_all if obs.area(th)[-1] >= 0.95 * obs.area(th).max() and th >= 0.1 and th < obs.ts[-1, 3]))

# ---------------- convergence table ----------------
conv = load_group("conv"); lines = [r"\begin{tabular}{cccccc}", r"\toprule", r"$L$ & $D_I$ & $A_I(T;0.01)$, $\Delta x = 0.04$ & $A_I(T;0.01)$, $\Delta x = 0.02$ & rel.\ diff.\ (\%) & regime \\", r"\midrule"]
maxpct = 0.0
for (L, D), rc in sorted(conv.items()):
    rf = fixed[(L, D)]; a_co, a_fi = rf.area(0.01)[-1], rc.area(0.01)[-1]; pct = 100 * abs(a_fi - a_co) / a_fi; maxpct = max(maxpct, pct)
    reg_c = classify(*rf.fractions(0.01, T_FIXED)); reg_f = classify(*rc.fractions(0.01, T_FIXED))
    lines.append(f"{L:g} & {D:g} & {a_co:.2f} & {a_fi:.2f} & {pct:.2f} & {reg_c}" + ("" if reg_c == reg_f else f" / {reg_f}") + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]; open(os.path.join(MS, "table_conv.tex"), "w").write("\n".join(lines) + "\n")
N("convmaxpct", f"{maxpct:.1f}"); print(f"convergence max rel diff {maxpct:.2f}%")

# ---------------- Fig 8: variants ----------------
local = load_group("local"); kernel = load_group("kernel"); T_LOC, T_KER = 25.0, 10.0
for d in sorted(glob.glob(os.path.join(DATA, "long", "*"))):
    if not os.path.exists(os.path.join(d, "summary.txt")): continue
    r = Run(d); a = r.area(0.01) / r.L ** 2; j = np.argmax(a >= 0.999) if np.any(a >= 0.999) else None; j5 = np.argmax(a >= 0.5)
    print(f"long {r.model}: 50% at t={r.t()[j5]:.0f}, 99.9% at t={r.t()[j] if j is not None else np.nan:.0f}, final {a[-1]:.3f}")
    N("longLocalT" if r.model == "local" else "longKernelT", f"{r.t()[j]:.0f}" if j is not None else "--")
    N("longLocalHalf" if r.model == "local" else "longKernelHalf", f"{r.t()[j5]:.0f}")
fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.8)); var_dc = {}; var_area = {}
for ax, name, runs, T in ((axs[0, 0], "local", local, T_LOC), (axs[0, 1], "kernel", kernel, T_KER)):
    Lv, Dv, labv, ffv, pkv = regime_table(runs, 0.01, T); dcv = crossover(Lv, Dv, ffv); var_dc[name] = dcv
    draw_map(ax, Lv, Dv, labv, title=f"({'ab'[name == 'kernel']}) {name} model, $T = {T:g}$", dots=dcv)
    Lf_, Df_ = inside(dcv); cntv = {r: int(np.sum(labv == r)) for r in REGIMES}; print(f"{name} T={T}: {cntv}; D_c inside grid:", dict(zip(Lf_, np.round(Df_, 4))), "censored:", [(L, dcv[L]) for L in dcv if dcv[L][1] != 0])
    key = "loc" if name == "local" else "ker"
    if len(Lf_) >= 2:
        p_, c_, r_ = powerfit(Lf_, Df_); N(key + "Exp", f"{p_:.1f}"); print(f"   power exponent {p_:.2f}")
        if len(Lf_) >= 3: co, l0, er = offsetfit(Lf_, Df_); print(f"   offset fit c={co:.2e} L0={l0:.2f} err {er:.1f}%")
    else: N(key + "Exp", "--")
    N(key + "Lmin", f"{Lf_.min():g}" if len(Lf_) else "--"); N(key + "Lmax", f"{Lf_.max():g}" if len(Lf_) else "--")
    i = len(Lv) - 1; m = (~np.isnan(pkv[i])) & (pkv[i] / Lv[i] ** 2 < 0.9); p_, c_, r_ = powerfit(Dv[m], pkv[i][m]); var_area[name] = (Dv[m], pkv[i][m], p_, c_, T)
    N(key + "Alpha", f"{p_:.2f}"); print(f"   area exponent at L={Lv[i]} (unsaturated points): {p_:.2f}")
    # intrinsic speed estimate at smallest D from the area at T
    a0 = pkv[i][0] / np.pi; N(key + "Speed", f"{np.sqrt(a0) / T:.2f}")
ax = axs[1, 0]
ax.plot(Lf, Df, "ko", ms=4, label="hotspot, $T = 400$"); ax.plot(Lline[Lline > L0 + 1], c_off * (Lline[Lline > L0 + 1] - L0) ** 2, "k-", lw=1.2)
for name, col, mk, T in (("local", "C0", "s", T_LOC), ("kernel", "C3", "^", T_KER)):
    dcv = var_dc[name]; Lf_, Df_ = inside(dcv)
    ax.plot(Lf_, Df_, mk, color=col, ms=5, label=f"{name}, $T = {T:g}$")
    Lc2 = [L for L in dcv if dcv[L][1] != 0]
    if Lc2: ax.plot(Lc2, [dcv[L][0] for L in Lc2], mk, mfc="w", mec=col, ms=5)
    if len(Lf_) >= 2: p_, c_, _ = powerfit(Lf_, Df_); ax.plot(Lline, c_ * Lline ** p_, "-", color=col, lw=1.2)
ax.plot(Lline, Lline ** 2 / (8 * np.pi * 1.0 * T_LOC ** 2), "--", color="C0", lw=1.2, label=rf"$L^2/(8\pi k_1 S_b T^2)$, $T = {T_LOC:g}$")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Domain size $L$"); ax.set_ylabel("Crossover boundary $D_c(L)$"); ax.set_title("(c) Crossover boundaries", fontsize=9)
ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both")
ax = axs[1, 1]; i = len(Ls) - 1; m = ~np.isnan(pk[i])
ax.plot(Ds[m], pk[i][m] / Ls[i] ** 2, "ko", ms=4, label=rf"hotspot, $T=400$: $D_I^{{{alphas[aL][0]:.2f}}}$"); xx = np.array([Ds[m][0], Ds[m][-1]]); ax.plot(xx, alphas[aL][1] * xx ** alphas[aL][0] / Ls[i] ** 2, "k-", lw=1)
for name, col, mk in (("local", "C0", "s"), ("kernel", "C3", "^")):
    Dv, y, p_, c_, T = var_area[name]; ax.plot(Dv, y / 18 ** 2, mk, color=col, ms=5, label=rf"{name}, $T={T:g}$: $D_I^{{{p_:.2f}}}$"); xx = np.array([Dv[0], Dv[-1]]); ax.plot(xx, c_ * xx ** p_ / 18 ** 2, "-", color=col, lw=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Diffusion coefficient $D_I$"); ax.set_ylabel(r"$A_I^{\rm peak}/L^2$ at $L = 18$"); ax.set_title("(d) Peak infected area fraction at $L = 18$", fontsize=9)
ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig9.eps")); fig.savefig(os.path.join(MS, "fig9.png")); plt.close(fig)

# ---------------- numbers.tex ----------------
for k in ("longLocalT", "longKernelT", "longLocalHalf", "longKernelHalf"): numbers.setdefault(k, "--")
with open(os.path.join(MS, "numbers.tex"), "w") as f:
    f.write("% data-dependent numbers, generated by code/analyze.py -- do not edit by hand\n")
    for k, v in numbers.items(): f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print("numbers:", numbers)
with open(os.path.join(DATA, "regimes_fixed_theta0.01.csv"), "w") as f:
    f.write("L,D,final_frac,peak_area,regime\n")
    for i, L in enumerate(Ls):
        for j, D in enumerate(Ds): f.write(f"{L:g},{D:g},{ff[i,j]:.4f},{pk[i,j]:.3f},{lab[i,j]}\n")
print("done")