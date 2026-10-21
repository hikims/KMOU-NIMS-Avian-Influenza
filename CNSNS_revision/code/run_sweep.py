#!/usr/bin/env python3
"""Run all simulations for the CNSNS revision (parallel driver).

Sweeps
  fixed  : localized-source model, corrected scheme, dx=0.04, dt=0.001, T=400
           L = 10..18, D = 0.001..0.010 (same grid as the submitted phase diagram)
  local  : local mass-action model (chi = 1),  dx=0.02, dt=0.0005, T=25
  kernel : nonlocal kernel model (r0=0.5, eps=0.1), dx=0.02, dt=0.0005, T=25
  conv   : fixed model at the submitted resolution dx=0.02, dt=0.0002 (convergence check)
  long   : local / kernel models at T=400 for one (L, D) (shows full saturation)

Usage:  python3 run_sweep.py [--workers 4] [--threads 8] [--only fixed,local,...]
"""
import argparse, itertools, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "sim")
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

ap = argparse.ArgumentParser()
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--only", type=str, default="")
ap.add_argument("--dry", action="store_true")
args = ap.parse_args()

jobs = []

def add(group, model, L, D, dx, dt, T, extra=()):
    out = os.path.join(DATA, group, f"L_{L:.1f}_D_{D:.4f}")
    cmd = [SIM, "--model", model, "--L", str(L), "--D", str(D), "--dx", str(dx), "--dt", str(dt),
           "--T", str(T), "--threads", str(args.threads), "--out", out] + list(extra)
    jobs.append((group, out, cmd))

# ---- A. fixed model, full sweep --------------------------------------------
L_fixed = [10, 11, 12, 13, 14, 15, 16, 17, 18]
D_fixed = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]
for L, D in itertools.product(L_fixed, D_fixed):
    add("fixed", "fixed", L, D, 0.04, 0.001, 400, ("--snap", "0,100,200,400"))

# ---- B/C. local and kernel models, short horizon -----------------------------
L_alt = [10, 12, 14, 16, 18]
D_alt = [0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1]
for L, D in itertools.product(L_alt, D_alt):
    snap = ("--snap", "0,5,10,25") if (L == 14 and D == 0.02) else ()
    add("local", "local", L, D, 0.02, 0.0005, 25, ("--out_every", "0.05") + snap)
    add("kernel", "kernel", L, D, 0.02, 0.0005, 25, ("--out_every", "0.05", "--eps", "0.1", "--r0", "0.5") + snap)

# ---- D. convergence check at the submitted resolution -------------------------
for L, D in [(10, 0.001), (15, 0.005), (18, 0.010)]:
    add("conv", "fixed", L, D, 0.02, 0.0002, 400, ("--snap", "0,100,200,400"))

# ---- E. long-horizon runs of the I-dependent models ---------------------------
add("long", "local", 18, 0.001, 0.02, 0.0005, 400, ("--out_every", "0.2"))
add("long", "kernel", 18, 0.001, 0.02, 0.0005, 400, ("--out_every", "0.2", "--eps", "0.1", "--r0", "0.5"))

if args.only:
    keep = set(args.only.split(","))
    jobs = [j for j in jobs if j[0] in keep]

# skip finished runs
def done(out):
    return os.path.exists(os.path.join(out, "summary.txt"))

todo = [j for j in jobs if not done(j[1])]
print(f"{len(jobs)} jobs, {len(todo)} to run, workers={args.workers}, threads={args.threads}", flush=True)
if args.dry:
    for g, out, cmd in todo: print(g, out)
    sys.exit(0)

# long jobs first so the tail of the schedule is short
order = {"conv": 0, "long": 1, "fixed": 2, "kernel": 3, "local": 4}
todo.sort(key=lambda j: order[j[0]])

def run(job):
    g, out, cmd = job
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(out, "log.txt"), "w") as lg:
        rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)
    print(f"[{g}] {os.path.basename(out)} rc={rc} {time.time()-t0:.0f}s", flush=True)
    return rc

t0 = time.time()
with ThreadPoolExecutor(max_workers=args.workers) as ex:
    rcs = list(ex.map(run, todo))
print(f"ALL DONE in {time.time()-t0:.0f}s, failures={sum(1 for r in rcs if r)}", flush=True)
