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
thresh = args.thresh

# calculate pairwise euclidean distances
# ingore cluster 11 cells in reference cell identity
not11 = reference_cell_ident[reference_cell_ident["Cluster"] != 11].index
is11 = reference_cell_ident[reference_cell_ident["Cluster"] == 11].index
distances = euclidean_distances(embeddings, embeddings.loc[not11])
distances = pd.DataFrame(
    distances,
    index=embeddings.index,
    columns=not11
)
# distances = euclidean_distances(embeddings, embeddings.loc[reference_cell_ident.index])
# distances = pd.DataFrame(
#     distances,
#     index=embeddings.index,
#     columns=reference_cell_ident.index
# )


# get the nearest reference cell for each cell
nearest_distances = distances.min(axis=1)
nearest_reference = distances.idxmin(axis=1)
thresh_distance = nearest_distances.quantile(thresh)


# assign cluster and color to each cell
cell_identity = pd.DataFrame(
    dict(
        Barcode=nearest_reference.index,
        Cluster=reference_cell_ident.loc[nearest_reference, "Cluster"].values,
        Color=reference_cell_ident.loc[nearest_reference, "Color"].values
    )
)
cell_identity.set_index("Barcode", inplace=True)
cell_identity.loc[is11, "Cluster"] = 11
cell_identity.loc[is11, "Color"] = "#D4D4D4"

# assign unique cluster and color to cells that are further than the threshold
unannotated = nearest_distances[nearest_distances > thresh_distance].index
cell_identity.loc[unannotated, "Cluster"] = 12
cell_identity.loc[unannotated, "Color"] = "#D4D4D4"
cell_identity.to_csv(args.cell_ident_csv, index=True)