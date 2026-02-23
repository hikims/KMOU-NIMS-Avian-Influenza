
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

PROJECT_PATH = os.environ.get("PROJECT_PATH", ".")
IMAGE_DIR = os.path.join(PROJECT_PATH, "image")
os.makedirs(IMAGE_DIR, exist_ok=True)

def xtick_formatter(x, pos):
    if abs(x) < 1e-10:
        return "0.0001"
    return f"{x:.3f}"

folder_list = [
    "di_0.0001","di_0.0025","di_0.0050","di_0.0075",
    "di_0.010","di_0.025","di_0.050","di_0.075",
    "di_0.100","di_0.125","di_0.150"
]

di_vals, peak_I = [], []

for folder in folder_list:
    fname = os.path.join(PROJECT_PATH, folder, "region.txt")
    if not os.path.exists(fname):
        continue
    df = pd.read_csv(fname, sep=r"\s+", header=None,
                     names=["t","I","S"], engine="python")
    di_vals.append(float(folder.split("_")[1]))
    peak_I.append(df["I"].max())

di_vals = np.array(di_vals)
peak_I = np.array(peak_I)
order = np.argsort(di_vals)

fig, ax = plt.subplots()
ax.plot(di_vals[order], peak_I[order],
        marker="o", linestyle="--", color="black")

ax.set_xlabel("Diffusion Coefficient")
ax.set_ylabel("Peak Infected Area")
ax.xaxis.set_major_formatter(mticker.FuncFormatter(xtick_formatter))

plt.tight_layout()
outname = os.path.join(IMAGE_DIR, "peak_I_vs_di.pdf")
fig.savefig(outname)
plt.close(fig)
