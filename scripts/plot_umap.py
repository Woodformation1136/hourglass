from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import sys

argv = sys.argv[1:]
if len(argv) < 3:
    print("Usage: plot_umap.py <umap_projection.csv> <clusters.csv> <output.png>")
    sys.exit(1)
    
# unpack arguments
umap_projection, clusters, output = argv

# load data
umap_projection = pd.read_csv(umap_projection)
clusters = pd.read_csv(clusters)
merged_data = pd.merge(umap_projection, clusters, on="Barcode")

# create subplots
fig, ax = plt.subplots(1, 1, figsize=(20, 15))

# plot umap projection
sns.scatterplot(
    x="UMAP-1",
    y="UMAP-2",
    hue="Cluster",
    data=merged_data,
    palette=[
        "#1F77B4",
        "#8C564B",
        "#FF7F0F",
        "#2AA02A",
        "#F8E71C",
        "#9467BD",
        "#D62728",
        "#E377C2",
        "#9B9B9B",
        "#4B4B4B"
    ],
    ax=ax
)

# save plot
plt.savefig(output)