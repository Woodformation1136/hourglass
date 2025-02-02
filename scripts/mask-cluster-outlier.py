from sklearn.metrics.pairwise import euclidean_distances
from pathlib import Path
import scipy.stats as stats
import pandas as pd
import numpy as np
import argparse
import os

# This script will mask outliers in clusters based on the correlation of intra-
# cluster distance on embedding 1 and 2

# create argument parser
parser = argparse.ArgumentParser(
    description="Get cell identity by nearest reference neighbor"
)
parser.add_argument(
    "--embedding1",
    type=str,
    required=True,
    help="The embeddings file where all cells are projected to (index format = sample:barcode)"
)
parser.add_argument(
    "--embedding2",
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
embedding1 = pd.read_csv(args.embedding1, index_col=0)
embedding2 = pd.read_csv(args.embedding2, index_col=0)
reference_cell_ident = pd.read_csv(args.reference_cell_ident, index_col=0)
out_cell_ident = reference_cell_ident.copy()
# reference_cell_ident["barcode"] = [i.split(":")[1] for i in reference_cell_ident.index]
thresh = args.thresh

# calculate pairwise euclidean distances
data1 = embedding1.merge(
    reference_cell_ident,
    left_index=True,
    right_index=True
)
data2 = embedding2.merge(
    reference_cell_ident,
    left_index=True,
    right_index=True
)
for cluster_n in reference_cell_ident["Cluster"].unique():
    # get the data for the cluster
    cluster_data1 = data1[data1["Cluster"] == cluster_n]
    cluster_data1 = cluster_data1.drop(columns=["Cluster", "Color"])
    cluster_data2 = data2[data2["Cluster"] == cluster_n]
    cluster_data2 = cluster_data2.drop(columns=["Cluster", "Color"])
    
    # calculate intra-cluster distances
    distances1 = pd.DataFrame(
        euclidean_distances(cluster_data1, cluster_data1),
        index=cluster_data1.index,
        columns=cluster_data1.index
    )
    distances2 = pd.DataFrame(
        euclidean_distances(cluster_data2, cluster_data2),
        index=cluster_data2.index,
        columns=cluster_data2.index
    )
    
    # calculate the correlation between distances1 and distance2 by row
    corr = []
    for i in distances1.index:
        corr.append(
            stats.pearsonr(
                distances1.loc[i],
                distances2.loc[i]
            ).statistic
        )
    corr = pd.DataFrame(
        corr,
        index=distances1.index,
        columns=["correlation"]
    )
    thresh_corr = corr["correlation"].quantile(1-thresh)
    out_cell_ident.loc[
        corr[corr["correlation"] < thresh_corr].index,
        "Cluster"
    ] = 11
    out_cell_ident.loc[
        corr[corr["correlation"] < thresh_corr].index,
        "Color"
    ] = "#000000"
    
out_cell_ident.to_csv(
    args.output_csv,
    index=True
)
