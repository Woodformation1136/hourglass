from matplotlib import pyplot as plt
import json
import seaborn as sns
import pandas as pd
import sys

argv = sys.argv[1:]
if len(argv) < 4:
    print("Usage: plot_umap.py <umap_projection.csv> <clusters.csv> <palette.json> <output.png>")
    sys.exit(1)
    
# unpack arguments
umap_projection, clusters, palette, output = argv

# load data
umap_projection = pd.read_csv(umap_projection)
clusters = pd.read_csv(clusters)
palette = json.load(open(palette))
palette = { int(k): v for k, v in palette.items() }
merged_data = pd.merge(umap_projection, clusters, on="Barcode")

# create subplots
fig, ax = plt.subplots(1, 1, figsize=(20, 15))

# plot umap projection
sns.scatterplot(
    x="UMAP-1",
    y="UMAP-2",
    hue="Cluster",
    data=merged_data,
    palette=palette,
    ax=ax
)

# save plot
plt.savefig(output)