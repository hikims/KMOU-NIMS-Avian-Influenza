
import os
import numpy as np
import matplotlib.pyplot as plt

PROJECT_PATH = os.environ.get("PROJECT_PATH", ".")
IMAGE_DIR = os.path.join(PROJECT_PATH, "paper_image")
os.makedirs(IMAGE_DIR, exist_ok=True)

filename = os.path.join(PROJECT_PATH, "data", "area_cluster.txt")
data = np.loadtxt(filename, comments="#")

t = data[:,1]
A = data[:,2:8].T

labels = [f"Cluster{i}" for i in range(6)]

plt.figure(figsize=(8,6))
plt.imshow(A, cmap="viridis", aspect="auto",
           origin="lower", extent=[t.min(), t.max(), 0, 5])

plt.colorbar(label="Cluster Area")
plt.yticks(range(6), labels)
plt.xlabel("Time")
plt.title("Cluster Area Heatmap")
plt.tight_layout()

outname = os.path.join(IMAGE_DIR, "cluster_area_heatmap.pdf")
plt.savefig(outname)
plt.close()
