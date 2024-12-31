from matplotlib import pyplot as plt
from pathlib import Path
import seaborn as sns
import pandas as pd
import sys
import os

argv = sys.argv[1:]
if len(argv) != 3:
    print("Usage: plot-umap.py <umap_projection.csv> <cell_identity.csv> <output_dir>")
    sys.exit(1)
    
# unpack arguments
umap_projection, cell_ident, output = argv

# load data
umap_projection = pd.read_csv(umap_projection)
cell_ident = pd.read_csv(cell_ident)
merged_data = pd.merge(umap_projection, cell_ident, on="Barcode")
merged_data["Sample"] = [x.split(":")[0] for x in merged_data["Barcode"]]

# create output directory
os.makedirs(output, exist_ok=True)

# create subplots
fig, ax = plt.subplots(1, 1, figsize=(20, 15))

# plot umap colored by Cluster
sns.scatterplot(
    x="UMAP-1",
    y="UMAP-2",
    hue="Cluster",
    data=merged_data.sample(frac=1, random_state=17),
    palette=dict(set(zip(merged_data["Cluster"], merged_data["Color"]))),
    s=10,
    ax=ax
)

# save plot
plt.savefig(Path(output) / "umap_hue=cluster.png")

# plot umap colored by sample
fig, ax = plt.subplots(1, 1, figsize=(20, 15))
sns.scatterplot(
    x="UMAP-1",
    y="UMAP-2",
    hue="Sample",
    data=merged_data.sample(frac=1, random_state=17),
    s=10,
    ax=ax
)
plt.savefig(Path(output) / "umap_hue=sample.png")