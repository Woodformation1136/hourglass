from matplotlib import pyplot as plt
from pathlib import Path
from argparse import ArgumentParser
import seaborn as sns
import pandas as pd
import json
import os

# create argument parser
argParser = ArgumentParser()
argParser.add_argument("--cell_embedding", help="Path to cell embedding CSV file")
argParser.add_argument("--cell_identity", help="Path to cell identity CSV file")
argParser.add_argument("--output_dir", help="Path to output directory", required=False)
argParser.add_argument("--palette", help="Path to palette JSON file (if provided override colors assigned to clusters)", required=False)
argParser.add_argument("--output_fpath", help="Output file name", required=False)
    
# parse arguments
args = argParser.parse_args()
cell_embedding = args.cell_embedding
cell_ident = args.cell_identity
output = args.output_dir
palette = args.palette
output_fpath = Path(args.output_fpath) if args.output_fpath else Path(output) / "umap_hue=cluster.png"

# load data
cell_embedding = pd.read_csv(cell_embedding)
cell_ident = pd.read_csv(cell_ident)
merged_data = pd.merge(cell_embedding, cell_ident, on="Barcode")
merged_data["Sample"] = [x.split(":")[0] for x in merged_data["Barcode"]]
merged_data["Cluster"] = merged_data["Cluster"].astype(str)
if palette:
    import json
    with open(palette, "r") as f:
        palette = json.load(f)

# create output directory
os.makedirs(output_fpath.parent, exist_ok=True)

# create subplots
fig, ax = plt.subplots(1, 1, figsize=(20, 15))

# plot umap colored by Cluster
sns.scatterplot(
    x="UMAP-1",
    y="UMAP-2",
    hue="Cluster",
    data=merged_data,
    palette=palette,
    s=50*merged_data.shape[0]/4000,
    ax=ax
)
plt.savefig(output_fpath.__str__())
