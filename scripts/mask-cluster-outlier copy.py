from sklearn.metrics.pairwise import euclidean_distances
from pathlib import Path
import pandas as pd
import numpy as np
import argparse
import os

# This script will mask outliers in clusters based on their distance to cluster median

# create argument parser
parser = argparse.ArgumentParser(
    description="Get cell identity by nearest reference neighbor"
)
parser.add_argument(
    "--embeddings",
    type=str,
    required=True,
    help="The embeddings file where all cells are projected to (index format = sample:barcode)"
)
parser.add_argument(
    "--reference_cell_ident",
    type=str,
    required=True,
    help="""\
    The cell identity of a subset of cells in the embedding file (index format = sample:barcode)
    """
)
parser.add_argument(
    "--output_csv",
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
embeddings = pd.read_csv(args.embeddings, index_col=0)
reference_cell_ident = pd.read_csv(args.reference_cell_ident, index_col=0)
out_cell_ident = reference_cell_ident.copy()
# reference_cell_ident["barcode"] = [i.split(":")[1] for i in reference_cell_ident.index]
thresh = args.thresh

# calculate pairwise euclidean distances
data = embeddings.merge(
    reference_cell_ident,
    left_index=True,
    right_index=True
)
for cluster_n in data["Cluster"].unique():
    cluster_data = data[data["Cluster"] == cluster_n]
    cluster_data = cluster_data.drop(columns=["Cluster", "Color"])
    cluster_median = cluster_data.median()
    dist2median = pd.DataFrame(
        euclidean_distances(cluster_data, [cluster_median]),
        index=cluster_data.index,
        columns=["distance"]
    )
    thresh_distance = dist2median["distance"].quantile(thresh)
    out_cell_ident.loc[
        dist2median[dist2median["distance"] > thresh_distance].index,
        "Cluster"
    ] = 11
    out_cell_ident.loc[
        dist2median[dist2median["distance"] > thresh_distance].index,
        "Color"
    ] = "#000000"
    
out_cell_ident.to_csv(
    args.output_csv,
    index=True
)
