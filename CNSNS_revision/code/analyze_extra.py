#!/usr/bin/env python3
"""Additional verification analyses for the revision (run after analyze.py):
  - extended convergence table (all runs in data/conv)
  - discrete mass-conservation error over all production runs
  - kernel refresh-period sensitivity (data/kerconv)
  - hotspot-size sensitivity of the crossover boundary and its offset L0 (data/hotsize)
  - horizon dependence of the boundary (hotspot: T = 200, 300, 400; local: T = 15, 20, 25)
  - coarse-vs-production comparison at high thresholds for (15, 0.01)
  - front-speed prediction with the half-area criterion for the local model
Writes fig9.eps/png, table_conv.tex (overwrites), and appends macros to numbers.tex.
"""
import glob, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import importlib.util, sys
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("an", os.path.join(HERE, "analyze.py"))
# import helpers only (analyze.py executes on import; we re-use its functions by exec of the helper part)
src = open(os.path.join(HERE, "analyze.py")).read().split("numbers = {}")[0]
exec(src)
DATA = os.path.normpath(os.path.join(HERE, "..", "data")); MS = os.path.normpath(os.path.join(HERE, "..", "manuscript"))
numbers = {}
def N(k, v): numbers[k] = v

fixed = load_group("fixed"); T_FIXED = 400.0
Ls, Ds, lab, ff, pk = regime_table(fixed, 0.01, T_FIXED); dc0 = crossover(Ls, Ds, ff); Lf0, Df0 = inside(dc0)
c0, L00, e0 = offsetfit(Lf0, Df0)

# ---------- 1. extended convergence table ----------
conv = load_group("conv")
lines = [r"\begin{tabular}{cccccc}", r"\toprule", r"$L$ & $D_I$ & $A_I(T;0.01)$, $\Delta x = 0.04$ & $A_I(T;0.01)$, $\Delta x = 0.02$ & rel.\ diff.\ (\%) & regime \\", r"\midrule"]
maxpct = 0.0; labels_same = True
for (L, D), rc in sorted(conv.items()):
    rf = fixed[(L, D)]; a_co, a_fi = rf.area(0.01)[-1], rc.area(0.01)[-1]; pct = 100 * abs(a_fi - a_co) / a_fi; maxpct = max(maxpct, pct)
    reg_c = classify(*rf.fractions(0.01, T_FIXED)); reg_f = classify(*rc.fractions(0.01, T_FIXED)); labels_same &= (reg_c == reg_f)
    pos = "invasion (far from boundary)" if (L, D) == (18.0, 0.001) else ("just below $D_c$" if (L, D) == (15.0, 0.003) else ("just above $D_c$" if (L, D) == (14.0, 0.003) else ""))
    lines.append(f"{L:g} & {D:g} & {a_co:.2f} & {a_fi:.2f} & {pct:.2f} & {reg_c}" + ("" if reg_c == reg_f else f" / {reg_f}") + (f" ({pos})" if pos else "") + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]; open(os.path.join(MS, "table_conv.tex"), "w").write("\n".join(lines) + "\n")
N("convmaxpct", f"{maxpct:.2f}"); N("convN", str(len(conv))); print(f"conv: {len(conv)} cases, max {maxpct:.3f}%, labels same: {labels_same}")
# boundary shift under refinement: interpolate D_c(15) from (15,0.003) fine & (15,0.004) coarse? only coarse grid available -> report fraction change instead
for (L, D), rc in conv.items():
    if (L, D) in ((15.0, 0.003), (14.0, 0.003)):
        print(f"  ({L:g},{D:g}) fraction coarse {fixed[(L,D)].fractions(0.01,T_FIXED)[1]:.4f} fine {rc.fractions(0.01,T_FIXED)[1]:.4f}")

# ---------- 2. mass conservation ----------
err = 0.0
for r in fixed.values():
    N0 = r.ts[0, 1] + r.ts[0, 2]; err = max(err, np.max(np.abs(r.ts[:, 1] + r.ts[:, 2] - N0) / N0))
N("massErr", f"{err:.1e}".replace("e-0", r"\times 10^{-").replace("e-", r"\times 10^{-") + "}"); print(f"max relative mass error over fixed runs: {err:.2e}")

# ---------- 3. kernel refresh period ----------
kc = {}
for d in sorted(glob.glob(os.path.join(DATA, "kerconv", "conv_every_*"))):
    if os.path.exists(os.path.join(d, "summary.txt")):
        r = Run(d); kc[int(d.split("_")[-1])] = r.area(0.01)[-1]
if kc:
    ref = kc[min(kc)]; dev = {k: 100 * abs(v - ref) / ref for k, v in kc.items()}; print("kernel refresh:", kc, dev)
    N("kerRefreshDev", f"{dev.get(10, np.nan):.2f}"); N("kerRefreshDevTwenty", f"{dev.get(20, np.nan):.2f}")

# ---------- 4. hotspot-size sensitivity ----------
variants = [("r0_0.25", 0.25, 1.0), ("eps_0.5", 0.5, 0.5), ("baseline", 0.5, 1.0), ("r0_1.0", 1.0, 1.0)]
hs = {}
coarse_base = load_group("fixed_hi_theta")
for tag, r0, eps in variants:
    runs = coarse_base if tag == "baseline" else load_group(os.path.join("hotsize", tag))
    if not runs: continue
    L_, D_, lab_, ff_, _ = regime_table(runs, 0.01, T_FIXED); dcv = crossover(L_, D_, ff_); Lf_, Df_ = inside(dcv)
    r05 = r0 + eps * np.arctanh(1 - 2 * 0.05)
    res = {"r0": r0, "eps": eps, "r05": r05, "Lf": Lf_, "Df": Df_, "nsat": int(np.sum(lab_ == "saturation"))}
    if len(Lf_) >= 3:
        c_, l0_, e_ = offsetfit(Lf_, Df_); p_, _, _ = powerfit(Lf_, Df_); res.update(c=c_, L0=l0_, err=e_, pexp=p_)
    hs[tag] = res; print(f"hotsize {tag}: r0={r0} eps={eps} r05={r05:.2f} nsat={res['nsat']} inside L={list(Lf_)}", {k: round(res[k], 3) for k in ("c", "L0", "err", "pexp") if k in res})
for tag, key in (("r0_0.25", "Small"), ("eps_0.5", "Sharp"), ("baseline", "Base"), ("r0_1.0", "Large")):
    if tag in hs and "L0" in hs[tag]:
        N("Lz" + key, f"{hs[tag]['L0']:.1f}"); N("Rf" + key, f"{hs[tag]['r05']:.1f}"); N("Pe" + key, f"{hs[tag]['pexp']:.1f}"); N("Ce" + key, fmt_e(hs[tag]["c"]))
    elif tag in hs:
        N("Lz" + key, "--"); N("Rf" + key, f"{hs[tag]['r05']:.1f}"); N("Pe" + key, "--"); N("Ce" + key, "--")

# ---------- 5. horizon dependence ----------
hor = {}
for T in (200.0, 300.0, 400.0):
    L_, D_, lab_, ff_, _ = regime_table(fixed, 0.01, T); dcv = crossover(L_, D_, ff_); Lf_, Df_ = inside(dcv)
    res = {"Lf": Lf_, "Df": Df_}
    if len(Lf_) >= 3: c_, l0_, e_ = offsetfit(Lf_, Df_); res.update(c=c_, L0=l0_, err=e_); print(f"hotspot T={T:g}: c={c_:.3e} L0={l0_:.2f} cT={c_*T:.3e} err {e_:.1f}% (L {Lf_.min():g}-{Lf_.max():g})")
    hor[T] = res
cT = [hor[T]["c"] * T for T in hor if "c" in hor[T]]
N("hotCTmin", fmt_e(min(cT))); N("hotCTmax", fmt_e(max(cT))); N("hotLzMin", f"{min(hor[T]['L0'] for T in hor if 'L0' in hor[T]):.1f}"); N("hotLzMax", f"{max(hor[T]['L0'] for T in hor if 'L0' in hor[T]):.1f}")
local = load_group("local"); horl = {}
for T in (15.0, 20.0, 25.0):
    L_, D_, lab_, ff_, _ = regime_table(local, 0.01, T); dcv = crossover(L_, D_, ff_); Lf_, Df_ = inside(dcv)
    res = {"Lf": Lf_, "Df": Df_}
    if len(Lf_) >= 2: p_, c_, _ = powerfit(Lf_, Df_); res.update(c=c_, p=p_); print(f"local T={T:g}: exp {p_:.2f} pref {c_:.3e} pref*T^2 {c_*T*T:.3e} (L {Lf_.min():g}-{Lf_.max():g})")
    horl[T] = res
cT2 = [horl[T]["c"] * T * T for T in horl if "c" in horl[T]]
N("locCTTmin", fmt_e(min(cT2))); N("locCTTmax", fmt_e(max(cT2)))
N("locPredHalfArea", fmt_e(1.0 / (8 * np.pi * 0.9 * 25.0 ** 2))); N("locPredWall", fmt_e(1.0 / (16 * 0.9 * 25.0 ** 2)))
# ratio measured / half-area prediction at T=25
if "c" in horl[25.0]:
    Lf_, Df_ = horl[25.0]["Lf"], horl[25.0]["Df"]; ratio = Df_ / (Lf_ ** 2 / (8 * np.pi * 0.9 * 625)); N("locRatioMin", f"{ratio.min():.2f}"); N("locRatioMax", f"{ratio.max():.2f}"); print("local measured / half-area prediction:", np.round(ratio, 2))

# ---------- 6. coarse vs production at high thresholds, (15, 0.01) ----------
obs = load_group("fixed_obs").get((15.0, 0.01)); co = coarse_base.get((15.0, 0.01))
if obs and co:
    rows = []
    for th in (0.3, 0.5, 0.6, 0.7, 0.8):
        lo = classify(*obs.fractions(th, T_FIXED)); lc = classify(*co.fractions(th, T_FIXED)); po, pc = obs.fractions(th, T_FIXED)[0], co.fractions(th, T_FIXED)[0]
        rows.append((th, lo, lc, po, pc)); print(f"  theta={th}: prod {lo} (peak {po:.3f})  coarse {lc} (peak {pc:.3f})")
    N("hiThetaAgree", "all five" if all(r[1] == r[2] for r in rows) else "not all"); N("hiThetaPeakDev", f"{max(abs(r[3]-r[4]) for r in rows):.3f}")

# ---------- Fig 9 ----------
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.5))
ax = axs[0]; Lline = np.linspace(9.5, 18.5, 60)
for tag, mk, col in (("r0_0.25", "v", "C0"), ("eps_0.5", "d", "C2"), ("baseline", "o", "k"), ("r0_1.0", "^", "C3")):
    if tag not in hs: continue
    h = hs[tag]; lab_ = rf"$r_0 = {h['r0']:g}$, $\epsilon = {h['eps']:g}$" + (rf" ($L_0 = {h['L0']:.1f}$)" if "L0" in h else "")
    ax.plot(h["Lf"], h["Df"], mk, color=col, ms=5, label=lab_)
    if "L0" in h: m = Lline > h["L0"] + 0.5; ax.plot(Lline[m], h["c"] * (Lline[m] - h["L0"]) ** 2, "-", color=col, lw=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.set_xlabel("Domain size $L$"); ax.set_ylabel("$D_c(L)$"); ax.set_title("(a) Hotspot size, $T = 400$", fontsize=9); ax.legend(fontsize=6.5); ax.grid(alpha=0.3, which="both")
ax = axs[1]
for T, mk in ((200.0, "s"), (300.0, "D"), (400.0, "o")):
    h = hor[T]; ax.plot(h["Lf"], h["Df"] * T, mk, color="k", ms=4, mfc="none" if T != 400 else "k", label=rf"hotspot, $T = {T:g}$")
for T, mk in ((15.0, "s"), (20.0, "D"), (25.0, "o")):
    h = horl[T]; ax.plot(h["Lf"], h["Df"] * T * T, mk, color="C0", ms=4, mfc="none" if T != 25 else "C0", label=rf"local, $T = {T:g}$")
ax.plot(Lline, Lline ** 2 / (8 * np.pi * 0.9), "--", color="C0", lw=1, label=r"$L^2/(8\pi k_1 S_0)$")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(Ls); ax.set_xticklabels([f"{L:g}" for L in Ls]); ax.set_xlabel("Domain size $L$"); ax.set_ylabel(r"$D_c\,T$ (hotspot),  $D_c\,T^2$ (local)")
ax.set_title("(b) Horizon dependence", fontsize=9); ax.legend(fontsize=6.5, ncol=2); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(MS, "fig4.eps")); fig.savefig(os.path.join(MS, "fig4.png")); plt.close(fig)

# ---------- append numbers ----------
existing = set(re.findall(r"\\newcommand\{\\(\w+)\}", open(os.path.join(MS, "numbers.tex")).read()))
with open(os.path.join(MS, "numbers.tex"), "a") as f:
    f.write("% --- analyze_extra.py ---\n")
    for k, v in numbers.items(): f.write(("\\renewcommand" if k in existing else "\\newcommand") + f"{{\\{k}}}{{{v}}}\n")
print("extra numbers:", numbers)
