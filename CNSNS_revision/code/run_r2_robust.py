#!/usr/bin/env python3
"""Robustness runs of the second revision (CNSNS-D-26-02308 R2), hotspot model, T = 400.

cons : the complete main sweep (L = 10..18, D = 0.001..0.010) at production resolution with the mass-conserving
       boundary rows (--consistent_bc); checks that no label and no boundary value depends on the 1e-3 boundary
       residual of the production scheme.
offc : off-centre source, (xc, yc) = (0.25, 0.25) and (0.25, 0.5), L = 12 and 16 (coarse resolution).
rect : rectangles (0, Lx) x (0, Lx/2), Lx = 12, 16, 20, 24, centred source (coarse resolution).
hoth : hotspot sizes r0 = 0.25 and 1.0 for the homogenization crossover, L = 12..28 (coarse resolution).
Coarse resolution: dx = 0.08, dt = 4e-3, mass-conserving rows (same as the extended sweep).
"""
import itertools, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = next((p for p in (os.path.join(HERE, "bin", "sim_new"), os.path.join(HERE, "bin", "sim"), os.path.join(HERE, "sim")) if os.path.exists(p)), os.path.join(HERE, "bin", "sim")); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TH = "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"
DG = [round(0.001 * 2 ** (k / 2), 6) for k in range(20)]
only = sys.argv[1].split(","); workers = int(sys.argv[2]); threads = int(sys.argv[3])
C = ["--model", "fixed", "--dx", "0.08", "--dt", "0.004", "--T", "400", "--thetas", TH, "--consistent_bc", "--threads", str(threads)]
jobs = []
if "offc" in only:
    for L, (xc, yc), D in itertools.product((12, 16), ((0.25, 0.25), (0.25, 0.5)), DG[:14]):
        jobs.append([SIM] + C + ["--L", str(L), "--D", str(D), "--xc", str(xc), "--yc", str(yc), "--out", os.path.join(DATA, "offc", f"xc{xc}_yc{yc}", f"L_{L:.1f}_D_{D:.6f}")])
if "rect" in only:
    for L, D in itertools.product((12, 16, 20, 24), DG):
        jobs.append([SIM] + C + ["--L", str(L), "--Ly", str(L / 2), "--D", str(D), "--out", os.path.join(DATA, "rect", f"L_{L:.1f}_D_{D:.6f}")])
if "hoth" in only:
    for r0, L, D in itertools.product((0.25, 1.0), (12, 16, 20, 24, 28), DG[5:]):
        jobs.append([SIM] + C + ["--L", str(L), "--D", str(D), "--r0", str(r0), "--out", os.path.join(DATA, "hoth", f"r0_{r0}", f"L_{L:.1f}_D_{D:.6f}")])
if "cons" in only:
    for L, D in itertools.product(range(10, 19), [round(0.001 * k, 3) for k in range(1, 11)]):
        jobs.append([SIM, "--model", "fixed", "--dx", "0.04", "--dt", "0.001", "--T", "400", "--thetas", TH, "--consistent_bc", "--threads", str(threads),
                     "--L", str(L), "--D", str(D), "--out", os.path.join(DATA, "cons", f"L_{L:.1f}_D_{D:.4f}")])
def run(cmd):
    out = cmd[-1]
    if os.path.exists(os.path.join(out, "summary.txt")): return 0
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg: rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"{out.split('/data/')[1]} rc={rc} {time.time()-t0:.0f}s", flush=True); return rc
print(f"{len(jobs)} jobs ({','.join(only)}), workers={workers}, threads={threads}", flush=True)
with ThreadPoolExecutor(max_workers=workers) as ex: rcs = list(ex.map(run, jobs))
print("ROBUST DONE failures=%d" % sum(1 for r in rcs if r), flush=True)
