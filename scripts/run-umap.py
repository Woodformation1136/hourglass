import pandas as pd
import json
import umap
import sys

argv = sys.argv[1:]
if len(argv) < 3:
    print("Usage: umap.py <pca_projections.csv> <umap_projections.csv> <params.json>  <random_state>")
    sys.exit(1)
    
pca_projections = pd.read_csv(argv[0], index_col=0)
params = json.load(open(argv[2]))
if len(argv) == 4:
    params.update({"random_state": int(argv[3])})
reducer = umap.UMAP(**params)
umap_projections = reducer.fit_transform(pca_projections)
out_df = pd.DataFrame(
    umap_projections,
    index=pca_projections.index,
    columns=["UMAP-1", "UMAP-2"]
)
out_df.index.name = "Barcode"
out_df.to_csv(argv[1])
