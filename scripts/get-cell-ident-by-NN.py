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
    "--cell_ident_csv",
    type=str,
    required=True,
    help="The output directory to save the results (split by sample)"
)

# parse arguments
args = parser.parse_args()
embeddings = pd.read_csv(args.embeddings, index_col=0)
reference_cell_ident = pd.read_csv(args.reference_cell_ident, index_col=0)

# calculate pairwise euclidean distances
distances = euclidean_distances(embeddings, embeddings)
distances = pd.DataFrame(
    distances,
    index=embeddings.index,
    columns=embeddings.index
)

# retain only the distances between the reference cells and all cells
distances = distances.loc[:, reference_cell_ident.index]

# get the nearest reference cell for each cell
nearest_reference = distances.idxmin(axis=1)

# assign cluster and color to each cell
cell_identity = pd.DataFrame(
    dict(
        Barcode=nearest_reference.index,
        Cluster=reference_cell_ident.loc[nearest_reference, "Cluster"].values,
        Color=reference_cell_ident.loc[nearest_reference, "Color"].values
    )
)
cell_identity.to_csv(args.cell_ident_csv, index=False)