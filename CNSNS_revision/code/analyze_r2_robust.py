#!/usr/bin/env python3
"""Robustness analyses of the second revision (run after analyze_r2.py).

 A. main sweep recomputed with mass-conserving boundary rows (data/cons) versus the production sweep
 B. coverage levels other than one half (0.3, 0.7): confined boundary versus free-space level curve
 C. horizon dependence of the homogenization crossover D_h (extended sweep, T' = 200, 300, 400)
 D. off-centre source (data/offc)
 E. rectangular domains (data/rect)
 F. hotspot size and the homogenization crossover (data/hoth)
 G. resolution of the radial free-space solver (data/free/fine)
Writes ../manuscript/numbers_r2_robust.tex, ../supplementary/S3_robustness.tex and figure robustness.eps/png.
"""
import glob, os, re, json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "analyze_r2.py")).read().split("# ==========================================================================================================\nfree = Free(\"base\")")[0]
exec(src)                                    # helpers: Run, load_group, Free, disk_in_square, classify, crossing, powerfit, offsetfit, ...
numbers = {}
def N(k, v): numbers[k] = v
free = Free("base"); fixed = load_group(os.path.join(R1, "fixed")); Ls, Ds, lab, ff, pk = regime_table(fixed, 0.01, T0)
rows_tex = []                                # (check, result) rows of Table S2

# ---------------- A. conservative rows ----------------
cons = load_group(os.path.join(DATA, "cons"), need_T=T0)
print("A. cons runs:", len(cons))
if len(cons) >= 90:
    dmax = 0.0; nlab = 0; dM = 0.0; drift = 0.0; ntot = 0; changes = []; dmain = 0.0
    for key, r in fixed.items():
        c = cons[key]; drift = max(drift, np.max(np.abs(c.ts[:, 1] + c.ts[:, 2] - c.ts[0, 1] - c.ts[0, 2])) / (c.ts[0, 1] + c.ts[0, 2]))
        dM = max(dM, abs(c.ts[-1, 1] / r.ts[-1, 1] - 1))
        for th in r.thetas:
            a, b = r.fractions(th, T0), c.fractions(th, T0); dmax = max(dmax, abs(a[1] - b[1])); nlab += classify(*a) != classify(*b); ntot += 1
            if th == 0.01: dmain = max(dmain, abs(a[1] - b[1]))
            if classify(*a) != classify(*b): changes.append((key, th, a[1], b[1], classify(*a), classify(*b)))
    Lc, Dc_, labc, ffc, pkc = regime_table(cons, 0.01, T0)
    dcp = {L: crossing(Ds, ff[i]) for i, L in enumerate(Ls)}; dcc = {L: crossing(Dc_, ffc[i]) for i, L in enumerate(Lc)}
    relb = max(abs(dcc[L][0] / dcp[L][0] - 1) for L in Ls if dcp[L][1] == 0)
    up3 = lambda v: f"{np.ceil(v * 1000 - 1e-9) / 1000:.3f}"     # upper bounds are rounded up
    N("consDa", up3(dmax)); N("consDaMain", up3(dmain)); N("consLab", nlab); print("   label changes:", changes); assert len(changes) == 1 and changes[0][1] != 0.01
    (kL, kD), kth, ka, kb, la, lb = changes[0]; N("consLabTheta", f"{kth:g}"); N("consLabL", f"{kL:g}"); N("consLabD", f"{kD:g}"); N("consLabA", f"{ka:.3f}"); N("consLabB", f"{kb:.3f}"); N("consLabFrom", la); N("consLabTo", lb); N("consN", ntot); N("consDM", up(100*dM, 1)); N("consDrift", fmt_e(drift, 0)); N("consBdy", up(100*relb, 3))
    print(f"   max |d a_T| {dmax:.4f} over {ntot} (run, threshold) pairs; label changes {nlab}; max dM_I {100*dM:.2f}%; drift {drift:.1e}; boundary change {100*relb:.3f}%")
    rows_tex.append(("Main sweep recomputed with mass-conserving boundary rows (90 runs, production resolution)",
                     f"mass error $\\le {fmt_e(drift, 0)}$; coverage changes by at most ${up3(dmax)}$ over {ntot} (run, threshold) pairs (${up3(dmain)}$ at $\\theta = 0.01$); $D_c(L)$ changes by at most ${up(100*relb, 3)}\\%$; $M_I(T)$ by at most ${up(100*dM, 1)}\\%$; one label changes, at $\\theta = {kth:g}$ in the run $(L, D_I) = ({kL:g}, {kD:g})$, whose coverage is ${ka:.3f}$ with the production rows ({la}) and ${kb:.3f}$ with the mass-conserving rows ({lb}); no label changes at $\\theta = 0.01$"))

else:
    for k in ("consDa", "consDaMain", "consLabTheta", "consLabL", "consLabD", "consLabA", "consLabB", "consLabFrom", "consLabTo", "consN", "consDM", "consDrift", "consBdy"): N(k, "{\\color{red}??}")      # placeholders until the sweep is complete

# ---------------- B. other coverage levels ----------------
res = {}
for level in (0.3, 0.5, 0.7):
    rel = []
    for i, L in enumerate(Ls):
        c = crossing(Ds, ff[i], level)
        if c[1] != 0: continue
        # free-space fraction of the square on the same grid
        fr = np.array([disk_in_square(free.A(D), L) / L ** 2 for D in Ds]); cf = crossing(Ds, fr, level)
        if cf[1] == 0: rel.append(100 * (c[0] / cf[0] - 1))
    res[level] = rel; print(f"B. level {level}: n={len(rel)}, rel diff range {min(rel):+.2f}..{max(rel):+.2f}%")
N("levLowRel", up(max(abs(v) for v in res[0.3]), 1)); N("levHighRel", up(max(abs(v) for v in res[0.7]), 1)); N("levHighSign", "lower" if np.mean(res[0.7]) < 0 else "higher")
rows_tex.append(("Coverage level $0.3$ and $0.7$ instead of $0.5$ (main sweep)", f"confined boundary within ${numbers['levLowRel']}\\%$ of the free-space level curve at level $0.3$; up to ${numbers['levHighRel']}\\%$ {numbers['levHighSign']} at level $0.7$, where the level set is close to the walls"))

# ---------------- C. horizon dependence of D_h ----------------
ext = load_group(os.path.join(DATA, "ext"), need_T=T0)
Le = sorted({k[0] for k in ext}); De = [round(0.001 * 2 ** (k / 2), 6) for k in range(20)]
def dh_at(runs, Lset, T, level=0.9, Dset=De):
    out = {}
    for L in Lset:
        d = [D for D in Dset if (L, D) in runs]
        if len(d) < 3: continue
        pr = np.array([runs[(L, D)].ts[runs[(L, D)].at(T), 5] for D in d]); c = crossing(d, pr, level)
        if c[1] == 0 and d[int(np.argmax(pr > level)) - 1] == Dset[Dset.index(d[int(np.argmax(pr > level))]) - 1]: out[L] = c[0]
    return out
from scipy.optimize import least_squares
HT = (150.0, 200.0, 250.0, 300.0, 350.0, 400.0)
def dh_points(level):
    pts = []
    for L in Le:
        d = [D for D in De if (L, D) in ext]
        if len(d) < 3: continue
        for T in HT:
            pr = np.array([ext[(L, D)].ts[ext[(L, D)].at(T), 5] for D in d]); c = crossing(d, pr, level); j = int(np.argmax(pr > level))
            if c[1] == 0 and j < len(d) - 1 and L <= 40: pts.append((L, T, c[0]))   # range of the fit of the manuscript; crossing not at the last simulated point
    return np.array(pts)
dfit = {}
for level in (0.8, 0.9, 0.95):
    pts = dh_points(level)
    f3 = lambda p: np.log(p[0] * (pts[:, 0] - p[1]) ** 2 / (np.pi ** 2 * (pts[:, 1] - p[2]))) - np.log(pts[:, 2])
    f2 = lambda p: np.log(p[0] * (pts[:, 0] - p[1]) ** 2 / (np.pi ** 2 * pts[:, 1])) - np.log(pts[:, 2])
    r3 = least_squares(f3, [0.8, 3.0, 50.0]); e3 = np.abs(np.exp(f3(r3.x)) - 1); r2 = least_squares(f2, [0.8, 3.0]); e2 = np.abs(np.exp(f2(r2.x)) - 1)
    dfit[level] = dict(pts=pts, p=r3.x, mean=100 * e3.mean(), max=100 * e3.max(), p2=r2.x, max2=100 * e2.max(), mean2=100 * e2.mean())
    print(f"C. PR={level}: {len(pts)} (L, T') points; D_h = {r3.x[0]:.2f} (L-{r3.x[1]:.2f})^2/(pi^2 (T'-{r3.x[2]:.0f})): mean err {100*e3.mean():.1f}%, max {100*e3.max():.0f}%; without delay: {r2.x[0]:.2f}, {r2.x[1]:.2f}, mean {100*e2.mean():.0f}%, max {100*e2.max():.0f}%")
d9 = dfit[0.9]; pts = d9["pts"]
N("DhDelayN", len(pts)); N("DhDelayC", f"{d9['p'][0]:.2f}"); N("DhDelayL", f"{d9['p'][1]:.1f}"); N("DhDelay", f"{d9['p'][2]:.0f}"); N("DhDelayMean", f"{d9['mean']:.0f}"); N("DhDelayMax", up(d9['max'], 0))
N("DhNoDelayMax", up(d9['max2'], 0)); N("DhDelayLow", f"{dfit[0.8]['p'][2]:.0f}"); N("DhDelayHigh", f"{dfit[0.95]['p'][2]:.0f}")
# product D_h T' at fixed L: change between T' = 200 and T' = 400
tab = {(L, T): D for L, T, D in pts}; drop = [100 * (tab[(L, 200.0)] * 200 / (tab[(L, 400.0)] * 400) - 1) for L in Le if (L, 200.0) in tab and (L, 400.0) in tab]
N("DhTdropMin", f"{min(drop):.0f}"); N("DhTdropMax", f"{max(drop):.0f}"); print(f"   D_h T' at T'=200 exceeds that at T'=400 by {min(drop):.0f}..{max(drop):.0f}% ({len(drop)} domain sizes)")
# time at which half of the final infected mass has been produced, in the runs next to D_h(L, 400)
half = []
for L in Le:
    if (L, 400.0) not in tab: continue
    d = [D for D in De if (L, D) in ext]; D = min(d, key=lambda x: abs(np.log(x / tab[(L, 400.0)]))); r = ext[(L, D)]; half.append(r.ts[int(np.argmax(r.ts[:, 1] >= 0.5 * r.ts[-1, 1])), 0])
N("DhHalfMin", f"{min(half):.0f}"); N("DhHalfMax", f"{max(half):.0f}"); print(f"   half-mass time in the runs next to D_h: {min(half):.0f}..{max(half):.0f}")
rows_tex.append(("Horizon dependence of the homogenization crossover ($T' = 150, \\dots, 400$, extended sweep)",
                 f"$D_h T'$ at $T' = 200$ exceeds its value at $T' = 400$ by ${numbers['DhTdropMin']}$ to ${numbers['DhTdropMax']}\\%$; with a delay, $D_h \\approx {numbers['DhDelayC']}\\,(L - {numbers['DhDelayL']})^2/[\\pi^2 (T' - {numbers['DhDelay']})]$ describes {len(pts)} values with a mean deviation of ${numbers['DhDelayMean']}\\%$ (maximum ${numbers['DhDelayMax']}\\%$); without the delay the maximum deviation is ${numbers['DhNoDelayMax']}\\%$"))

# ---------------- C2. free-space control of the homogenization crossover (numbers from analyze_r2.py) ----------------
_m = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{([^{}]*)\}", open(os.path.join(MS, "numbers_r2.tex")).read()))
rows_tex.append(("Free-space control for the homogenization crossover: participation ratio of the free-space solution restricted to the square (extended grid)",
                 f"reaches $0.9$ for $L \\le {_m['DhFreeLmax']}$ at ${_m['DhFreeRatioMin']}$ to ${_m['DhFreeRatioMax']}$ times the confined $D_h$ ($L \\le 24$; ${_m['DhFreeRatioLast']}$ times at $L = {_m['DhFreeLmax']}$), and not within the sampled range $D_I \\le 0.72$ for $L \\ge {_m['DhFreeLnone']}$; factors ${_m['DhFreeRatioLowMin']}$ to ${_m['DhFreeRatioLowMax']}$ at level $0.8$ and ${_m['DhFreeRatioHighMin']}$ to ${_m['DhFreeRatioHighMax']}$ at level $0.95$; largest difference of participation ratio ${_m['prWallMax']}$"))

# ---------------- D. off-centre source ----------------
def frac_disk_offcentre(A, L, xc, yc, n=2400):
    r = np.sqrt(max(A, 0.0) / np.pi); x = (np.arange(n) + 0.5) * L / n; X, Y = np.meshgrid(x, x, indexing="ij")
    return float(np.mean((X - xc * L) ** 2 + (Y - yc * L) ** 2 <= r * r))
off = []
for tag, (xc, yc) in (("xc0.25_yc0.25", (0.25, 0.25)), ("xc0.25_yc0.5", (0.25, 0.5))):
    runs = load_group(os.path.join(DATA, "offc", tag), need_T=T0)
    if len(runs) < 28: continue
    for L in (12.0, 16.0):
        d = sorted(D for (l, D) in runs if l == L); cov = np.array([runs[(L, D)].area(0.01)[-1] / L ** 2 for D in d])
        covf = np.array([frac_disk_offcentre(free.A(D), L, xc, yc) for D in d]); c = crossing(d, cov); cf = crossing(d, covf); ccen = free.Dc(L)
        j = int(np.argmax(cov > 0.5)); wall = max(abs(cov[j] - covf[j]), abs(cov[j - 1] - covf[j - 1])); mrel = max(abs(runs[(L, D)].ts[-1, 1] / free.M(D) - 1) for D in d[: j + 1])
        drift = max(np.max(np.abs(r.ts[:, 1] + r.ts[:, 2] - r.ts[0, 1] - r.ts[0, 2])) / (r.ts[0, 1] + r.ts[0, 2]) for r in runs.values())
        print(f"D. off-centre {tag} L={L:g}: D_c={c[0]:.5f} (free-space disk clipped by the square: {cf[0]:.5f}, centred source: {ccen:.5f}); wall effect at boundary {100*wall:.1f}% of domain; mass diff {100*mrel:.1f}%; drift {drift:.0e}")
        off.append(dict(ratio=c[0] / ccen, rel=100 * (c[0] / cf[0] - 1), wall=100 * wall, mass=100 * mrel))
def s1(q, x0): return 1 + 2 * sum(np.cos(k * np.pi * x0) ** 2 * np.exp(-2 * k * k * q) for k in range(1, 200))     # point release at x0 L: 1/PR = s1(q, x0) s1(q, y0), q = D pi^2 T / L^2
from scipy.optimize import brentq as _brentq
qc = _brentq(lambda q: 1 / (s1(q, 0.5) ** 2) - 0.9, 1e-3, 50); offdh = {}
for tag, (xc, yc), key in (("xc0.25_yc0.25", (0.25, 0.25), "Corner"), ("xc0.25_yc0.5", (0.25, 0.5), "Side")):
    runs = load_group(os.path.join(DATA, "offc", tag), need_T=T0)
    if len(runs) < 28: continue
    d = sorted(D for (l, D) in runs if l == 12.0); c = crossing(d, np.array([runs[(12.0, D)].ts[-1, 5] for D in d]), 0.9); cen = dh_at(ext, [12.0], T0)[12.0]
    qo = _brentq(lambda q: 1 / (s1(q, xc) * s1(q, yc)) - 0.9, 1e-3, 80); assert c[1] == 0
    _r, _I = free.prof[De[0]]; _n = 600; _x = (np.arange(_n) + 0.5) * 12.0 / _n; _X, _Y = np.meshgrid(_x, _x, indexing="ij"); _R = np.hypot(_X - xc * 12.0, _Y - yc * 12.0)
    def _prf(D):                                                 # participation ratio of the free-space solution (source at the same position) restricted to the square
        r_, I_ = free.prof[D]; F_ = np.interp(_R, r_, I_, right=0.0); return F_.sum() ** 2 / (_n * _n * (F_ * F_).sum())
    cfr = crossing(De, np.array([_prf(D) for D in De]), 0.9); assert cfr[1] == 0; N("offFree" + key, f"{cfr[0] / c[0]:.1f}"); print(f"   free-space solution restricted to the square reaches PR = 0.9 at {cfr[0]:.4f} = {cfr[0]/c[0]:.2f} x confined")
    offdh[key] = (c[0] / cen, qo / qc); N("offDh" + key, f"{c[0] / cen:.1f}"); N("offPt" + key, f"{qo / qc:.1f}")
    print(f"D. off-centre {tag}, L=12: D_h = {c[0]:.4f} = {c[0]/cen:.2f} x centred; point-release estimate {qo/qc:.2f} (centred q = {qc:.3f})")
if len(off) == 4:
    N("offRatioMin", f"{min(o['ratio'] for o in off):.1f}"); N("offRatioMax", f"{max(o['ratio'] for o in off):.1f}"); N("offRelMin", f"{min(abs(o['rel']) for o in off):.1f}"); N("offRelMax", f"{max(abs(o['rel']) for o in off):.1f}")
    N("offWallMin", f"{min(o['wall'] for o in off):.1f}"); N("offWallMax", f"{max(o['wall'] for o in off):.1f}"); N("offMass", f"{max(o['mass'] for o in off):.0f}"); assert all(o["rel"] < 0 for o in off)
    rows_tex.append(("Off-centre source at $(L/4, L/4)$ and $(L/4, L/2)$, $L = 12, 16$ (56 runs)",
                     f"$D_c$ is ${numbers['offRatioMin']}$ to ${numbers['offRatioMax']}$ times that of the centred source, because part of the free-space disk falls outside the square; the confined boundary lies ${numbers['offRelMin']}$ to ${numbers['offRelMax']}\\%$ below the diffusivity at which the free-space solution, clipped by the square, covers half of it, and the wall effect at the boundary is ${numbers['offWallMin']}$ to ${numbers['offWallMax']}\\%$ of the domain (centred source: at most $0.3\\%$); at $L = 12$ the homogenization crossover is ${numbers['offDhCorner']}$ and ${numbers['offDhSide']}$ times that of the centred source (point-release estimate: ${numbers['offPtCorner']}$ and ${numbers['offPtSide']}$), and the free-space solution restricted to the square reaches the same participation ratio at ${numbers['offFreeCorner']}$ and ${numbers['offFreeSide']}$ times these diffusivities (centred source: about $2$)"))

# ---------------- E. rectangles ----------------
def offsetfit_rel(L, D):
    f = lambda p: np.log(p[0] * (L - p[1]) ** 2 / (np.pi ** 2 * T0)) - np.log(D); r = least_squares(f, [0.8, 3.0]); return r.x[0], r.x[1], 100 * np.max(np.abs(np.exp(f(r.x)) - 1))
def frac_disk_rect(A, Lx, Ly, n=3200):
    r = np.sqrt(max(A, 0.0) / np.pi); x = (np.arange(n) + 0.5) * Lx / n; y = (np.arange(n // 2) + 0.5) * Ly / (n // 2); X, Y = np.meshgrid(x, y, indexing="ij")
    return float(np.mean((X - Lx / 2) ** 2 + (Y - Ly / 2) ** 2 <= r * r))
def prq(qx, qy): return 1.0 / ((1 + 2 * sum(np.exp(-8 * qx * j * j) for j in range(1, 40))) * (1 + 2 * sum(np.exp(-8 * qy * j * j) for j in range(1, 40))))
from scipy.optimize import brentq
q_sq = brentq(lambda q: prq(q, q) - 0.9, 1e-3, 5); q_re = brentq(lambda q: prq(q, 4 * q) - 0.9, 1e-3, 5)     # point release: square of side L_x, rectangle L_x x L_x/2
rect = load_group(os.path.join(DATA, "rect"), need_T=T0)
print("E. rect runs:", len(rect), f"| point-release estimate of D_h(rect)/D_h(square of side L_x) = {q_re/q_sq:.2f}")
if len(rect) >= 80:
    out = []; drift = max(np.max(np.abs(r.ts[:, 1] + r.ts[:, 2] - r.ts[0, 1] - r.ts[0, 2])) / (r.ts[0, 1] + r.ts[0, 2]) for r in rect.values())
    for L in (12.0, 16.0, 20.0, 24.0):
        d = sorted(D for (l, D) in rect if l == L); pr = np.array([rect[(L, D)].ts[-1, 5] for D in d]); ch = crossing(d, pr, 0.9)
        cov = np.array([rect[(L, D)].area(0.01)[-1] / (L * L / 2) for D in d]); cc = crossing(d, cov); covf = np.array([frac_disk_rect(free.A(D), L, L / 2) for D in d]); cf = crossing(d, covf)
        j = int(np.argmax(cov > 0.5)); wall = max(abs(cov[j] - covf[j]), abs(cov[j - 1] - covf[j - 1])) if cc[1] == 0 else np.nan
        sq = dh_at(ext, [L], T0).get(L, np.nan); sqc = crossing(De, [ext[(L, D)].area(0.01)[-1] / L ** 2 if (L, D) in ext else disk_in_square(free.A(D), L) / L ** 2 for D in De])[0]
        out.append((L, ch[0] if ch[1] == 0 else np.nan, sq, cc[0] if cc[1] == 0 else np.nan, cf[0] if cf[1] == 0 else np.nan, wall, sqc))
        print(f"   Lx={L:g}, Ly={L/2:g}: D_h={out[-1][1]:.4f} (square of side Lx: {sq:.4f}, ratio {out[-1][1]/sq:.2f}); D_c={out[-1][3]:.5f} (free-space disk clipped by the rectangle: {out[-1][4]:.5f}, {100*(out[-1][3]/out[-1][4]-1):+.1f}%; square of side Lx: {sqc:.5f}); wall effect at boundary {100*wall:.1f}%")
    out = np.array(out); ratio = out[:, 1] / out[:, 2]; rel = 100 * (out[:, 3] / out[:, 4] - 1); print(f"   drift {drift:.0e}")
    N("rectRatioMin", f"{np.nanmin(ratio):.2f}"); N("rectRatioMax", f"{np.nanmax(ratio):.2f}"); N("rectPt", f"{q_re/q_sq:.2f}"); N("rectRelMin", f"{np.nanmin(np.abs(rel)):.1f}"); N("rectRelMax", f"{np.nanmax(np.abs(rel)):.1f}"); assert np.nanmax(rel) < 0
    sh = [out[k, 1] / dh_at(ext, [out[k, 0] / 2], T0)[out[k, 0] / 2] for k in range(len(out)) if out[k, 0] / 2 in Le]; N("rectShortMin", f"{min(sh):.0f}"); N("rectShortMax", f"{max(sh):.0f}"); print("   D_h(rect) / D_h(square of the shorter side):", np.round(sh, 2))
    N("rectWallMin", f"{np.nanmin(out[:, 5])*100:.1f}"); N("rectWallMax", f"{np.nanmax(out[:, 5])*100:.1f}"); N("rectDcRatioMin", f"{np.nanmin(out[:, 3]/out[:, 6]):.2f}"); N("rectDcRatioMax", f"{np.nanmax(out[:, 3]/out[:, 6]):.2f}")
    rows_tex.append(("Rectangles $L_x \\times L_x/2$, $L_x = 12, 16, 20, 24$, centred source (80 runs)",
                     f"$D_h$ is ${numbers['rectRatioMin']}$ to ${numbers['rectRatioMax']}$ times its value in the square of side $L_x$ (point-release estimate: ${numbers['rectPt']}$) and ${numbers['rectShortMin']}$ to ${numbers['rectShortMax']}$ times its value in the square of side $L_x/2$ ($L_x = 20, 24$), i.e.\\ it is governed by the longer side; the coverage boundary lies ${numbers['rectRelMin']}$ to ${numbers['rectRelMax']}\\%$ below the diffusivity at which the free-space solution, clipped by the rectangle, covers half of it, with a wall effect of ${numbers['rectWallMin']}$ to ${numbers['rectWallMax']}\\%$ of the domain at the boundary"))

# ---------------- F. hotspot size and D_h ----------------
hoth = {}; LF = [12.0, 16.0, 20.0, 24.0, 28.0]
for r0 in (0.25, 1.0):
    runs = load_group(os.path.join(DATA, "hoth", f"r0_{r0}"), need_T=T0)
    if len(runs) < 75: continue
    dh = dh_at(runs, LF, T0, Dset=De[5:]); Lh = np.array(sorted(dh)); Dh = np.array([dh[L] for L in Lh])
    if len(Lh) >= 4: c_, l0_, e_ = offsetfit_rel(Lh, Dh); hoth[r0] = (c_, l0_, e_, Lh, Dh); print(f"F. r0={r0}: D_h {dict(zip(Lh, np.round(Dh, 4)))}; fit {c_:.2f} (L-{l0_:.2f})^2/(pi^2T), max err {e_:.1f}%")
dhb = dh_at(ext, LF, T0); Lb = np.array(sorted(dhb)); Db = np.array([dhb[L] for L in Lb]); cb, l0b, eb = offsetfit_rel(Lb, Db)
print(f"   baseline r0=0.5 on the same L: D_h {dict(zip(Lb, np.round(Db, 4)))}; fit {cb:.2f} (L-{l0b:.2f})^2/(pi^2T), max err {eb:.1f}%")
if len(hoth) == 2:
    N("LhSmall", f"{hoth[0.25][1]:.1f}"); N("LhBase", f"{l0b:.1f}"); N("LhLarge", f"{hoth[1.0][1]:.1f}"); N("ChSmall", f"{hoth[0.25][0]:.2f}"); N("ChBase", f"{cb:.2f}"); N("ChLarge", f"{hoth[1.0][0]:.2f}")
    rs = [hoth[0.25][4][list(hoth[0.25][3]).index(L)] / Db[list(Lb).index(L)] for L in Lb if L in hoth[0.25][3]]; rl = [hoth[1.0][4][list(hoth[1.0][3]).index(L)] / Db[list(Lb).index(L)] for L in Lb if L in hoth[1.0][3]]
    N("DhSmallRatioMin", f"{min(rs):.2f}"); N("DhSmallRatioMax", f"{max(rs):.2f}"); N("DhLargeRatioMin", f"{min(rl):.2f}"); N("DhLargeRatioMax", f"{max(rl):.2f}"); print("   D_h(r0)/D_h(0.5): small", np.round(rs, 3), "large", np.round(rl, 3))
    rows_tex.append(("Hotspot radius $r_0 = 0.25$, $0.5$, $1$ and the homogenization crossover ($L = 12, \\dots, 28$; 150 runs)",
                     f"$D_h$ changes by a factor of ${numbers['DhSmallRatioMin']}$ to ${numbers['DhSmallRatioMax']}$ ($r_0 = 0.25$) and ${numbers['DhLargeRatioMin']}$ to ${numbers['DhLargeRatioMax']}$ ($r_0 = 1$) relative to the baseline; fitted offsets $L_h = {numbers['LhSmall']}$, ${numbers['LhBase']}$, ${numbers['LhLarge']}$ and prefactors ${numbers['ChSmall']}$, ${numbers['ChBase']}$, ${numbers['ChLarge']}$"))

# ---------------- H. diffusivities beyond the extended sweep (data/extD): the crossing guaranteed as D -> infinity ----------------
xd = load_group(os.path.join(DATA, "extD"), need_T=T0); LD = (44.0, 48.0, 56.0)
if len(xd) >= 18:
    prmin, dmin, dcross, mlarge, mlast, dr_ = [], [], [], [], [], []
    for L in LD:
        runs = {D: r for (l, D), r in list(ext.items()) + list(xd.items()) if l == L}; d = np.array(sorted(runs)); pr = np.array([runs[D].ts[-1, 5] for D in d]); M = np.array([runs[D].ts[-1, 1] for D in d])
        k0 = int(np.argmax(pr[d <= De[-1]])); k1 = k0 + int(np.argmin(pr[k0:])); up_ = [k for k in range(k1, len(d) - 1) if pr[k] < 0.9 <= pr[k + 1]]; k = up_[0]
        dc_ = float(np.exp(np.log(d[k]) + (0.9 - pr[k]) * (np.log(d[k + 1]) - np.log(d[k])) / (pr[k + 1] - pr[k])))
        prmin.append(pr[k1]); dmin.append(d[k1]); dcross.append(dc_); mlarge.append(M[d > De[-1]].max()); mlast.append(M[d == De[-1]][0]); assert np.all(np.array([runs[D].area(0.01)[-1] for D in d if D > De[-1]]) == 0)
        dr_.append(max(np.max(np.abs(r.ts[:, 1] + r.ts[:, 2] - r.ts[0, 1] - r.ts[0, 2])) / (r.ts[0, 1] + r.ts[0, 2]) for (l, D), r in xd.items() if l == L))
        print(f"H. L={L:g}: PR minimum {pr[k1]:.3f} at D={d[k1]:.3f}, high-D upward crossing of 0.9 at D={dc_:.3f}; max M_I(T) for D > {De[-1]:.2f}: {mlarge[-1]:.2f} (at D = {De[-1]:.2f}: {mlast[-1]:.2f}); PR at largest D {pr[-1]:.4f}")
    lst = lambda v, f: ", ".join(f.format(x) for x in v[:-1]) + " and " + f.format(v[-1])
    lstm = lambda v, f: "$, $".join(f.format(x) for x in v[:-1]) + "$ and $" + f.format(v[-1])      # the same inside $...$
    N("DhLDn", len(xd)); N("DhLDdrift", fmt_e(max(dr_), 0)); N("DhLDmass", up(max(mlarge), 1)); N("DhLDmassSweep", f"{min(mlast):.1f}")
    # resolution checks of the two runs bracketing the crossing at L = 44 (data/extDchk): halved time step; halved grid spacing with a quarter time step
    chk = {tag: load_group(os.path.join(DATA, "extDchk", tag), need_T=T0) for tag in ("dt_half", "dx_half")}; chknote = ""; chkrow = ""
    if all(len(v) == 2 for v in chk.values()):
        base = {D: r for (l, D), r in xd.items() if l == 44.0}; dd = sorted(D for (l, D) in chk["dt_half"]); res = {}
        for tag, runs in chk.items():
            dpr = [abs(runs[(44.0, D)].ts[-1, 5] - base[D].ts[-1, 5]) for D in dd]; dm = [abs(runs[(44.0, D)].ts[-1, 1] / base[D].ts[-1, 1] - 1) for D in dd]
            prc = [runs[(44.0, D)].ts[-1, 5] for D in dd]; dc2 = float(np.exp(np.log(dd[0]) + (0.9 - prc[0]) * (np.log(dd[1]) - np.log(dd[0])) / (prc[1] - prc[0])))
            res[tag] = (max(dpr), 100 * max(dm), 100 * abs(dc2 / dcross[0] - 1)); print(f"   check {tag}: max |dPR| {max(dpr):.5f}, max rel dM {100*max(dm):.3f}%, crossing {dc2:.4f} ({100*(dc2/dcross[0]-1):+.2f}%)")
        N("DhLDchkPR", up(max(v[0] for v in res.values()), 3)); N("DhLDchkM", up(max(v[1] for v in res.values()), 1)); N("DhLDchkDc", up(max(v[2] for v in res.values()), 1))
        chknote = f" Repeating the two runs that bracket the crossing at $L = 44$ with half the time step, and with half the grid spacing and a quarter of the time step, changes $\\mathrm{{PR}}(T)$ by at most ${numbers['DhLDchkPR']}$, the infected mass by at most ${numbers['DhLDchkM']}\\%$ (the coarse grid resolves the seed with about three cells per radius) and the interpolated crossing by at most ${numbers['DhLDchkDc']}\\%$."
        chkrow = f"; with half the time step, and with half the grid spacing ($L = 44$, two diffusivities), $\\mathrm{{PR}}(T)$ changes by at most ${numbers['DhLDchkPR']}$, $M_I(T)$ by at most ${numbers['DhLDchkM']}\\%$ and the crossing by at most ${numbers['DhLDchkDc']}\\%$"
    N("DhLargeDNote", f"Additional runs for $L = 44$, $48$ and $56$ at diffusivities up to ${max(D for (l, D) in xd):.1f}$ (Supplementary Material) show this high-diffusivity branch: beyond the sweep the participation ratio falls to a minimum of ${lstm(prmin, '{:.2f}')}$ near $D_I \\approx {np.mean(dmin):.1f}$ and then rises towards one, crossing $0.9$ at $D_I \\approx {lstm(dcross, '{:.1f}')}$ (values interpolated on the grid of diffusivities). At these diffusivities much less infected mass is produced than inside the sweep ($M_I(T) \\le {up(max(mlarge), 1)}$, against at least ${min(mlast):.1f}$ at $D_I = 0.72$), the coverage at the horizon is zero, and the crossing depends only weakly on $L$. This return to uniformity at large diffusivity lies outside the range described by \\eqref{{eq:Dh}}, and we have not analysed it further.{chknote}")
    N("DhLargeDLetter", f"In {len(xd)} additional runs for L = 44, 48 and 56 at diffusivities up to {max(D for (l, D) in xd):.1f} the participation ratio falls to a minimum near a diffusivity of {np.mean(dmin):.1f} and then rises towards one, crossing 0.9 at about {lst(dcross, '{:.1f}')}; much less infected mass is produced there than inside the sweep and the coverage at the horizon is zero. This return to uniformity at large diffusivity lies outside the range described by the fit, and we have not analysed it further.")
    rows_tex.append((f"Diffusivities beyond the extended sweep ($D_I$ up to ${max(D for (l, D) in xd):.1f}$; $L = 44, 48, 56$; {len(xd)} runs)",
                     f"$\\mathrm{{PR}}(T)$ has a minimum of ${lstm(prmin, '{:.2f}')}$ near $D_I \\approx {np.mean(dmin):.1f}$ and then tends to one, crossing $0.9$ at $D_I \\approx {lstm(dcross, '{:.2f}')}$; there $M_I(T) \\le {up(max(mlarge), 1)}$ and the coverage at the horizon is zero; mass error $\\le {fmt_e(max(dr_), 0)}${chkrow}"))
else:
    N("DhLargeDNote", ""); N("DhLargeDLetter", "")

# ---------------- G. radial solver resolution ----------------
fine = Free("fine")
if len(fine.Dg) >= 3:
    dA = max(abs(fine.ts[D][-1, 7] / free.ts[D][-1, 7] - 1) for D in fine.Dg if D in free.ts); dM = max(abs(fine.ts[D][-1, 1] / free.ts[D][-1, 1] - 1) for D in fine.Dg if D in free.ts)
    N("radA", up(100*dA, 3)); N("radM", up(100*dM, 3)); N("radN", len(fine.Dg)); print(f"G. radial solver, dr and dt halved ({len(fine.Dg)} values of D): max dA {100*dA:.3f}%, dM {100*dM:.3f}%")
    def _prsq(prof, L, n=600):
        r_, I_ = prof; x_ = (np.arange(n) + 0.5) * L / n - L / 2; X_, Y_ = np.meshgrid(x_, x_, indexing="ij"); F_ = np.interp(np.hypot(X_, Y_), r_, I_, right=0.0); return F_.sum() ** 2 / (n * n * (F_ * F_).sum())
    dprp = max(abs(_prsq(fine.prof[D], L) - _prsq(free.prof[D], L)) for D in fine.Dg if D in free.prof for L in (10.0, 16.0, 24.0, 40.0)); N("radPR", up(dprp, 4))
    print(f"   free-space PR in the square, refined versus base radial solution (stored spacing 0.05 vs 0.10): max |dPR| {dprp:.2e}")
    rows_tex.append((f"Radial solver with $\\Delta r$ and $\\Delta t$ halved ({len(fine.Dg)} diffusivities from $10^{{-3}}$ to $0.13$)", f"thresholded area changes by at most ${up(100*dA, 3)}\\%$, infected mass by at most ${up(100*dM, 3)}\\%$, participation ratio of the solution restricted to the square by at most ${up(dprp, 4)}$"))

with open(os.path.join(MS, "numbers_r2_robust.tex"), "w") as f:
    f.write("% robustness numbers of the second revision, generated by code/analyze_r2_robust.py\n")
    for k, v in numbers.items(): f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
with open(os.path.join(SUP, "S3_robustness.tex"), "w") as f:
    f.write("\\begin{tabular}{p{6.2cm}p{9.2cm}}\n\\toprule\nCheck & Result \\\\\n\\midrule\n")
    for a, b in rows_tex: f.write(f"{a} & {b} \\\\\n\\addlinespace\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
json.dump(rows_tex, open(os.path.join(SUP, "S3_robustness.json"), "w"), indent=1)
print("robust numbers:", numbers)
