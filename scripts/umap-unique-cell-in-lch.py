from sklearn.metrics.pairwise import euclidean_distances
from pathlib import Path
import pandas as pd
import argparse
import os

# create argument parser
parser = argparse.ArgumentParser(
    description="Get cell identity by nearest reference neighbor"
)
parser.add_argument(
    "--umap_embedding",
    type=str,
    required=True,
    help="The embeddings file where all cells are projected to (index format = sample:barcode)"
)
parser.add_argument(
    "--ref_cells",
    type=str,
    required=True
)
parser.add_argument(
    "--input_cell_ident",
    type=str,
    required=True
)
parser.add_argument(
    "--output_cell_ident",
    type=str,
    required=True,
    help="The output directory to save the results (split by sample)"
)
parser.add_argument(
    "--thresh",
    type=float,
    default=1,
    help="The NN further then the qth-percentile of all distance will be annotated as unique cells"
)

# parse arguments
args = parser.parse_args()
embeddings = pd.read_csv(args.umap_embedding, index_col=0)
ref_cells = list(open(args.ref_cells, "r").read().splitlines())
input_cell_ident = pd.read_csv(args.input_cell_ident, index_col=0)
thresh = args.thresh

# calculate pairwise euclidean distances
# ingore cluster 11 cells in reference cell identity
not11_all_cells = input_cell_ident[input_cell_ident["Cluster"] != 11].index
is11_all_cells = input_cell_ident[input_cell_ident["Cluster"] == 11].index
not11_ref_cells = list(set(ref_cells) & set(not11_all_cells))
not_ref_cells = list(set(embeddings.index) - set(ref_cells))
distances = euclidean_distances(embeddings.loc[not_ref_cells], embeddings.loc[not11_ref_cells])
distances = pd.DataFrame(
    distances,
    index=embeddings.loc[not_ref_cells].index,
    columns=not11_ref_cells
)


# identity unique cells that are further than the threshold
nearest_distances = distances.min(axis=1)
nearest_reference = distances.idxmin(axis=1)
thresh_distance = nearest_distances.quantile(thresh)
unique_cells = nearest_distances[nearest_distances > thresh_distance].index
lch_green = list(
    set(input_cell_ident[input_cell_ident["Cluster"] == 5].index) &
    set(input_cell_ident[input_cell_ident.index.map(lambda x: "lch" in x)].index)
)
lch_blue = list(
    set(input_cell_ident[input_cell_ident["Cluster"] == 1].index) &
    set(input_cell_ident[input_cell_ident.index.map(lambda x: "lch" in x)].index)
)
lch_brown = list(
    set(input_cell_ident[input_cell_ident["Cluster"] == 2].index) &
    set(input_cell_ident[input_cell_ident.index.map(lambda x: "lch" in x)].index)
)
lch_fusiform_late = list(set(lch_brown) | set(lch_blue))

unique_lch_green = list(set(unique_cells) & set(lch_green))
unique_lch_fusiform_late = list(set(unique_cells) & set(lch_fusiform_late))

# update input_cell_identity
input_cell_ident.loc[unique_lch_green, "Cluster"] = 13
input_cell_ident.loc[unique_lch_green, "Color"] = "#33C7FF"
input_cell_ident.loc[unique_lch_fusiform_late, "Cluster"] = 14
input_cell_ident.loc[unique_lch_fusiform_late, "Color"] = "#33C7FF"
input_cell_ident.to_csv(args.output_cell_ident, index=True)
