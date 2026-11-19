#!/usr/bin/env python3
"""Second-revision runs (CNSNS-D-26-02308 R2).

ext   : extended (L, D_I) sweep of the hotspot model, L = 10..56, D_I = 0.001..0.724 (factor sqrt(2) steps),
        coarse resolution dx = 0.08, dt = 4e-3, T = 400, mass-conserving boundary rows (--consistent_bc).
        Purpose: sample domain sizes well beyond the fitted offset L0 (reviewer point 2) and the regime in
        which the walls matter (reviewer point 4).
ref2d : two-dimensional reference runs in a large domain (L = 40) at production resolution, used to validate
        the radial free-space solver radial_free.py.
extchk: production-resolution checks (dx = 0.04, dt = 1e-3) of extended-sweep runs next to the boundary.
"""
import itertools, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = next((p for p in (os.path.join(HERE, "bin", "sim_new"), os.path.join(HERE, "bin", "sim"), os.path.join(HERE, "sim")) if os.path.exists(p)), os.path.join(HERE, "bin", "sim")); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TH = "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"
L_EXT = [10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56]
D_EXT = [round(0.001 * 2 ** (k / 2), 6) for k in range(20)]          # 0.001 ... 0.724
only = set(sys.argv[1].split(",")) if len(sys.argv) > 1 else {"ext"}
workers = int(sys.argv[2]) if len(sys.argv) > 2 else 5
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
jobs = []
if "ext" in only:
    # Runs in which the cloud stays far from the walls coincide with the free-space reference (radial_free.py)
    # and are skipped: a pair (L, D) is run only if the free-space density at the distance L/2 of the nearest
    # wall exceeds 1e-5 at t = T (list written by the selection step into ext_needed.json).
    import json
    needed = sorted([tuple(p) for p in json.load(open(os.path.join(HERE, "ext_needed.json")))], key=lambda p: (p[0], p[1]))
    for L, D in needed:
        out = os.path.join(DATA, "ext", f"L_{L:.1f}_D_{D:.6f}")
        jobs.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.08", "--dt", "0.004", "--T", "400",
                     "--threads", str(threads), "--thetas", TH, "--consistent_bc", "--out", out])
if "ref2d" in only:
    for D in (0.005, 0.01):
        out = os.path.join(DATA, "ref2d", f"L_40.0_D_{D:.6f}")
        jobs.append([SIM, "--model", "fixed", "--L", "40", "--D", str(D), "--dx", "0.04", "--dt", "0.001", "--T", "400",
                     "--threads", str(threads), "--thetas", TH, "--out", out])
if "extchk" in only:
    for L, D in [(20, 0.008), (20, 0.011314), (24, 0.016), (24, 0.022627), (20, 0.064)]:   # next to D_c (first four) and to D_h (last)
        out = os.path.join(DATA, "extchk", f"L_{L:.1f}_D_{D:.6f}")
        jobs.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.04", "--dt", "0.001", "--T", "400",
                     "--threads", str(threads), "--thetas", TH, "--consistent_bc", "--out", out])
def run(cmd):
    out = cmd[-1]
    if os.path.exists(os.path.join(out, "summary.txt")): return 0
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg: rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"{out.split('/data/')[1]} rc={rc} {time.time()-t0:.0f}s", flush=True); return rc
print(f"{len(jobs)} jobs ({','.join(sorted(only))}), workers={workers}, threads={threads}", flush=True)
with ThreadPoolExecutor(max_workers=workers) as ex: rcs = list(ex.map(run, jobs))
print("R2 RUNS DONE failures=%d" % sum(1 for r in rcs if r), flush=True)
