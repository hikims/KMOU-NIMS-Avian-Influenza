#!/usr/bin/env python3
"""Boundary-treatment and kernel-refresh checks reported in Section 4.1 of the revised manuscript.

(a) Boundary rows: six (L, D_I) pairs at the coarse resolution (dx = 0.08, dt = 4e-3, T = 400) run twice,
    with the mixed boundary rows of the production scheme and with --consistent_bc ((-1,1) rows in both
    sweeps, mass conserved to 1e-12). The set contains the wall-contact case (10, 0.01), an invasion case far
    from the boundary and the grid points on both sides of the interpolated crossover boundary at L = 14, 15.
(b) Kernel refresh interval: the representative kernel run (14, 0.02, T = 25) and two runs on opposite
    sides of the kernel boundary at T = 10 (L = 14: D = 0.05 invasion, D = 0.07 saturation), each with the
    convolution refreshed every step and every ten steps.

Outputs: ../data/bc_check/{mixed,consistent}/L_*_D_*/ and ../data/kerconv_T10/conv_every_{1,10}/L_14.0_D_*/
Compare final areas / masses with:  python3 run_bc_check.py --report
"""
import os, subprocess, sys, itertools
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.join(HERE, "sim"); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
PAIRS = [(10, 0.01), (18, 0.001), (15, 0.003), (15, 0.004), (14, 0.002), (14, 0.003)]
TH = "0.005,0.01,0.02,0.05,0.1,0.2,0.3"

def jobs():
    out = []
    for L, D in PAIRS:
        for tag, extra in (("mixed", []), ("consistent", ["--consistent_bc"])):
            d = os.path.join(DATA, "bc_check", tag, f"L_{L:.1f}_D_{D:.4f}")
            out.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.08", "--dt", "0.004", "--T", "400", "--threads", "4", "--thetas", TH, "--out", d] + extra)
    for D in (0.02, 0.05, 0.07):
        for ce in (1, 10):
            d = os.path.join(DATA, "kerconv_T10", f"conv_every_{ce}", f"L_14.0_D_{D:.4f}")
            out.append([SIM, "--model", "kernel", "--L", "14", "--D", str(D), "--dx", "0.02", "--dt", "0.0005", "--T", "25", "--threads", "4", "--eps", "0.1", "--r0", "0.5", "--out_every", "0.05", "--conv_every", str(ce), "--out", d])
    return out

def run(cmd):
    d = cmd[-1]
    if os.path.exists(os.path.join(d, "summary.txt")): return 0
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "log.txt"), "w") as lg: return subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)

def report():
    import numpy as np
    def final(d, col):
        a = np.loadtxt(os.path.join(d, "timeseries.txt"), comments="#"); return a[-1, col], a
    print("boundary rows: L D | area mixed / consistent (rel %) | massI mixed / consistent (rel %) | mass drift mixed / consistent")
    for L, D in PAIRS:
        dm = os.path.join(DATA, "bc_check", "mixed", f"L_{L:.1f}_D_{D:.4f}"); dc = os.path.join(DATA, "bc_check", "consistent", f"L_{L:.1f}_D_{D:.4f}")
        if not (os.path.exists(dm + "/summary.txt") and os.path.exists(dc + "/summary.txt")): continue
        am, A = final(dm, 8); ac, B = final(dc, 8); mm, _ = final(dm, 1); mc, _ = final(dc, 1)
        drift = lambda a: np.max(np.abs(a[:, 1] + a[:, 2] - (a[0, 1] + a[0, 2])) / (a[0, 1] + a[0, 2]))
        print(f"  {L:2d} {D:.3f} | {am:8.2f} {ac:8.2f} ({100*abs(am-ac)/ac:.2f}%) | {mm:7.2f} {mc:7.2f} ({100*abs(mm-mc)/mc:.2f}%) | {drift(A):.1e} {drift(B):.1e}")
    print("kernel refresh at T=10 (coverage at t=10) and T=25:")
    for D in (0.02, 0.05, 0.07):
        res = {}
        for ce in (1, 10):
            d = os.path.join(DATA, "kerconv_T10", f"conv_every_{ce}", f"L_14.0_D_{D:.4f}")
            if os.path.exists(d + "/summary.txt"):
                a = np.loadtxt(d + "/timeseries.txt", comments="#"); j = np.argmin(np.abs(a[:, 0] - 10)); res[ce] = (a[j, 8] / 196, a[-1, 8] / 196)
        if res: print(f"  D={D}: " + ", ".join(f"every {ce}: t=10 {v[0]:.4f}, t=25 {v[1]:.4f}" for ce, v in res.items()))

if __name__ == "__main__":
    if "--report" in sys.argv: report(); sys.exit(0)
    with ThreadPoolExecutor(max_workers=4) as ex: rcs = list(ex.map(run, jobs()))
    print("done, failures:", sum(1 for r in rcs if r)); report()
