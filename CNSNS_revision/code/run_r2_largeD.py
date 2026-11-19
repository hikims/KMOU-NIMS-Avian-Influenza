#!/usr/bin/env python3
"""Second revision, check of the large-diffusivity limit: runs beyond the diffusivity range of the extended sweep
(D = 1.024 ... 5.793 in steps of sqrt(2)) for L = 44, 48, 56, to locate the homogenization crossover PR(T) = 0.9 on the
high-diffusivity branch (PR(T) -> 1 as D -> infinity by estimate (19)). Same grid as the extended sweep (dx = 0.08, mass-conserving boundary rows);
the time step is reduced so that D dt / dx^2 <= 1.   Output: ../data/extD/L_*_D_*/   Usage: python3 run_r2_largeD.py [workers] [threads] [--checks]   (--checks: only the four resolution checks, ../data/extDchk/)"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
SIM = next((p for p in (os.path.join(HERE, "bin", "sim_new"), os.path.join(HERE, "bin", "sim"), os.path.join(HERE, "sim")) if os.path.exists(p)), os.path.join(HERE, "bin", "sim"))
W = int(sys.argv[1]) if len(sys.argv) > 1 else 6; TH = sys.argv[2] if len(sys.argv) > 2 else "4"
jobs = []
for L in (44, 48, 56):
    for k in range(20, 26):
        D = round(0.001 * 2 ** (k / 2), 6); dt = 0.004 if D <= 1.5 else (0.002 if D <= 3.0 else 0.001); assert D * dt / 0.08 ** 2 <= 1.0
        jobs.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.08", "--dt", str(dt), "--T", "400", "--consistent_bc", "--threads", TH,
                     "--thetas", "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8", "--out", os.path.join(DATA, "extD", f"L_{L:.1f}_D_{D:.6f}")])
def run(cmd):
    out = cmd[-1]
    if os.path.exists(os.path.join(out, "summary.txt")): return 0
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg: rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"{out.split('/data/')[1]} rc={rc} {time.time()-t0:.0f}s", flush=True); return rc
# resolution checks of the pair of runs bracketing the high-diffusivity crossing at L = 44: halved time step, and halved grid spacing with a quarter time step
if "--checks" in sys.argv:
    jobs = []
    for D in (1.024, 1.448155):
        for tag, dx, dt, th in (("dt_half", "0.08", "0.002", TH), ("dx_half", "0.04", "0.001", "12")):
            jobs.append([SIM, "--model", "fixed", "--L", "44", "--D", str(D), "--dx", dx, "--dt", dt, "--T", "400", "--consistent_bc", "--threads", th,
                         "--thetas", "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8", "--out", os.path.join(DATA, "extDchk", tag, f"L_44.0_D_{D:.6f}")])
print(f"{len(jobs)} jobs (extD), workers={W}, threads={TH}", flush=True)
with ThreadPoolExecutor(W) as ex: rcs = list(ex.map(run, sorted(jobs, key=lambda c: (float(c[6]), float(c[4])))))
print("LARGE-D DONE failures=%d" % sum(r != 0 for r in rcs), flush=True)
