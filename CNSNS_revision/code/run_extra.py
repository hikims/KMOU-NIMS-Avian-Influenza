#!/usr/bin/env python3
"""Extra runs for the revision: (a) hotspot sweep at coarse resolution with extended thresholds
(for the high-threshold regime maps), (b) one representative run with extended thresholds at
production resolution (observable comparison)."""
import itertools, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.join(HERE, "sim"); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TH = "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"
jobs = []
out = os.path.join(DATA, "fixed_obs", "L_15.0_D_0.0100")
jobs.append([SIM, "--model", "fixed", "--L", "15", "--D", "0.01", "--dx", "0.04", "--dt", "0.001", "--T", "400", "--threads", "4", "--thetas", TH, "--out", out])
for L, D in itertools.product(range(10, 19), [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]):
    out = os.path.join(DATA, "fixed_hi_theta", f"L_{L:.1f}_D_{D:.4f}")
    jobs.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.08", "--dt", "0.004", "--T", "400", "--threads", "4", "--thetas", TH, "--out", out])
def run(cmd):
    out = cmd[-1]; os.makedirs(out, exist_ok=True); t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg: rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"{os.path.basename(os.path.dirname(out))}/{os.path.basename(out)} rc={rc} {time.time()-t0:.0f}s", flush=True)
with ThreadPoolExecutor(max_workers=4) as ex: list(ex.map(run, jobs))
print("EXTRA DONE", flush=True)
