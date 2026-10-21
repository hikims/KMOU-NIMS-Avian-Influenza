#!/usr/bin/env python3
"""Second batch of extra runs (internal review, 2026-10-15):
  conv2    : fine-resolution convergence checks on invasion and near-boundary cases
  hotsize  : hotspot-size sensitivity (r0 = 0.25, 1.0 and eps = 0.5) at coarse resolution, full (L, D) grid
  kerconv  : kernel refresh-period sensitivity (conv_every = 1, 5, 10)
"""
import itertools, os, subprocess, time, sys
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.join(HERE, "sim"); DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TH = "0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"
slow, fast = [], []
# conv2: invasion case, just below and just above the boundary
for L, D in [(18, 0.001), (15, 0.003), (14, 0.003)]:
    out = os.path.join(DATA, "conv", f"L_{L:.1f}_D_{D:.4f}")
    slow.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.02", "--dt", "0.0002", "--T", "400", "--threads", "8", "--out", out])
# hotspot-size sensitivity (coarse grid validated against production: identical labels)
for tag, r0, eps in (("r0_0.25", 0.25, 1.0), ("r0_1.0", 1.0, 1.0), ("eps_0.5", 0.5, 0.5)):
    for L, D in itertools.product(range(10, 19), [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]):
        out = os.path.join(DATA, "hotsize", tag, f"L_{L:.1f}_D_{D:.4f}")
        fast.append([SIM, "--model", "fixed", "--L", str(L), "--D", str(D), "--dx", "0.08", "--dt", "0.004", "--T", "400", "--threads", "4", "--r0", str(r0), "--eps", str(eps), "--thetas", TH, "--out", out])
# kernel refresh-period sensitivity
for ce in (1, 5, 10, 20):
    out = os.path.join(DATA, "kerconv", f"conv_every_{ce}")
    fast.append([SIM, "--model", "kernel", "--L", "14", "--D", "0.02", "--dx", "0.02", "--dt", "0.0005", "--T", "25", "--threads", "4", "--eps", "0.1", "--r0", "0.5", "--out_every", "0.05", "--conv_every", str(ce), "--out", out])
def run(cmd):
    out = cmd[-1]
    if os.path.exists(os.path.join(out, "summary.txt")): return 0
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg: rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"{out.split('/data/')[1]} rc={rc} {time.time()-t0:.0f}s", flush=True); return rc
with ThreadPoolExecutor(max_workers=5) as ex:
    fs = [ex.submit(run, c) for c in slow] + [ex.submit(run, c) for c in fast]
    rcs = [f.result() for f in fs]
print("EXTRA2 DONE failures=%d" % sum(1 for r in rcs if r), flush=True)
