from argparse import ArgumentParser
import pandas as pd
import json
import umap
import os

# create argument parser
argParser = ArgumentParser()
argParser.add_argument("--projections_in", help="Path to PCA projections CSV file")
argParser.add_argument("--projections_out", help="Path to output UMAP projections CSV file")
argParser.add_argument("--params", help="Path to UMAP parameters JSON file")
argParser.add_argument("--random_state", help="Random state for UMAP")
argParser.add_argument("--n_components", help="Number of components for UMAP", required=False)
argParser.add_argument("--cell_identity", help="Cell identity file name", required=False)
argParser.add_argument("--target_clusters", help="Target cluster for UMAP", required=False)

# parse arguments
args = argParser.parse_args()
projections_in = args.projections_in
projections_out = args.projections_out
params = args.params
random_state = args.random_state
n_components = args.n_components
cell_identity = args.cell_identity
target_clusters = args.target_clusters

# load data    
projections_in = pd.read_csv(projections_in, index_col=0)
params = json.load(open(params))
if random_state:
    params.update({"random_state": int(random_state)})
if n_components:
    params.update({"n_components": int(n_components)})
if target_clusters:
    if cell_identity is None:
        raise ValueError("Cell identity file name is required to specify target cluster")
    cell_identity = pd.read_csv(cell_identity)
    target_clusters = set([int(i) for i in target_clusters.split(",")])
    target_cells = cell_identity[cell_identity["Cluster"].isin(target_clusters)]["Barcode"].values
    projections_in = projections_in.loc[target_cells]

    
# run UMAP
print("Running UMAP with parameters:")
print(params)
reducer = umap.UMAP(**params)
umap_projections = reducer.fit_transform(projections_in)
out_df = pd.DataFrame(
    umap_projections,
    index=projections_in.index,
    columns=[f"UMAP-{i+1}" for i in range(int(n_components))]
)
out_df.index.name = "Barcode"
os.makedirs(os.path.dirname(projections_out), exist_ok=True)
out_df.to_csv(projections_out)
