
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

PROJECT_PATH = os.environ.get("PROJECT_PATH", ".")
IMAGE_DIR = os.path.join(PROJECT_PATH, "image")
os.makedirs(IMAGE_DIR, exist_ok=True)

folder_list = [
    "di_0.0001","di_0.0025","di_0.0050","di_0.0075",
    "di_0.010","di_0.025","di_0.050","di_0.075",
    "di_0.100","di_0.125","di_0.150"
]

for folder in folder_list:
    fname = os.path.join(PROJECT_PATH, folder, "region.txt")
    if not os.path.exists(fname):
        print(f"[WARNING] {fname} not found.")
        continue

    df = pd.read_csv(fname, sep=r"\s+", header=None,
                     names=["t","I","S"], engine="python")

    df_plot = df.iloc[::100, :].copy()

    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.plot(df_plot["t"], df_plot["I"], color="red", label="Infected Area")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Infected Area", color="red")

    ax2 = ax1.twinx()
    ax2.plot(df_plot["t"], df_plot["S"], color="blue",
             linestyle="--", label="Susceptible Area")
    ax2.set_ylabel("Susceptible Area", color="blue")

    plt.title("Time Series of Infected and Susceptible Areas")
    plt.tight_layout()

    outname = os.path.join(IMAGE_DIR, f"{folder}_time_series_IS.pdf")
    fig.savefig(outname)
    plt.close(fig)

print("Completed.")
