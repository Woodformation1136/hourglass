from matplotlib import pyplot as plt
from pathlib import Path
from argparse import ArgumentParser
import seaborn as sns
import pandas as pd
import numpy as np
import json
import os

# create argument parser
argParser = ArgumentParser()
argParser.add_argument("--cell_embedding", help="Path to cell embedding CSV file")
argParser.add_argument("--palette", help="Path to palette JSON file (if provided override colors assigned to clusters)", required=False)
argParser.add_argument("--output_prefix", help="Output file name", required=False)
argParser.add_argument("--n_components", help="Number of components for UMAP", required=False, default=2)
    
# parse arguments
args = argParser.parse_args()
cell_embedding = args.cell_embedding
palette = args.palette
output_prefix = Path(args.output_prefix)
n_components = int(args.n_components)

# load data
cell_embedding = pd.read_csv(cell_embedding)
cell_embedding["Sample"] = [i.split(":")[0] for i in cell_embedding["Barcode"]]
if palette:
    import json
    with open(palette, "r") as f:
        palette = json.load(f)

# create output directory
os.makedirs(output_prefix.parent, exist_ok=True)

# plot facet grid of umap colored by Cluster
# row = "Species"
# col = "Batch"
if n_components == 2:
    # plot hue = sample
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(
        data=cell_embedding,
        x="UMAP-1",
        y="UMAP-2",
        hue="Sample",
        s=1,
        ax=ax
    )
    fig.savefig(
        output_prefix.with_suffix("").__str__() + "_hue=sample.png",
        dpi=300
    )
    # # plot hue = cluster
    # g = sns.FacetGrid(
    #     merged_data,
    #     row=row,
    #     col=col,
    #     hue="Cluster",
    #     palette=palette,
    #     height=5,
    #     aspect=1
    # )
    # g.map_dataframe(
    #     sns.scatterplot,
    #     x="UMAP-1",
    #     y="UMAP-2",
    #     s=10
    # )
    # plt.savefig(
    #     output_prefix.with_suffix("").__str__() + "_hue=cluster.png",
    #     dpi=300
    # )

# # test elevation and azimuth using "ptr" data
# if n_components == 3:
#     from mpl_toolkits.mplot3d import axes3d
#     elevations = np.linspace(0, 300, 6)
#     azimuths = np.linspace(0, 300, 6)
    
#     # configure plot
#     plt.style.use("fast")
#     fig, axes = plt.subplots(
#         nrows=len(elevations),
#         ncols=len(azimuths),
#         figsize=(40,40),
#         subplot_kw=dict(projection="3d")
#     )
#     for nth_row in range(len(elevations)):
#         for nth_col in range(len(azimuths)):
#             # plot hue = sample
#             ax = axes[nth_row, nth_col]
#             r_name = elevations[nth_row]
#             c_name = azimuths[nth_col]
#             subset = merged_data[merged_data["Sample"] == "ptr_tenx_batch1_rs17"]
#             ax.view_init(elev=r_name, azim=c_name, roll=0)
#             # ax.view_init(elev=60, azim=-30, roll=0) # TODO: try this view
#             if subset.shape[0] == 0:
#                 continue
#             ax.scatter(
#                 subset["UMAP-1"],
#                 subset["UMAP-2"],
#                 subset["UMAP-3"],
#                 color=subset["Cluster"].map(palette),
#                 s=1
#             )
#             ax.set_xlabel("UMAP-1")
#             ax.set_ylabel("UMAP-2")
#             ax.set_zlabel("UMAP-3")
#             ax.set_title(f"elevation={r_name}, azimuths={c_name}")
#             ax.xaxis.pane.fill = False
#             ax.yaxis.pane.fill = False
#             ax.zaxis.pane.fill = False

#             # set color to white (or whatever is "invisible")
#             ax.xaxis.pane.set_edgecolor('w')
#             ax.yaxis.pane.set_edgecolor('w')
#             ax.zaxis.pane.set_edgecolor('w')

#             # get rid of the grid as well:
#             ax.grid(linestyle=":")
            
#     # save plot
#     plt.savefig(
#         output_prefix.with_suffix("").__str__() + "_hue=cluster.png",
#         dpi=200
#     )
    
# # plot by Species
# elev, azim, roll = 120, 180, 0
# if n_components == 3:
#     from mpl_toolkits.mplot3d import axes3d
#     row_names = merged_data[row].unique()
#     col_names = merged_data[col].unique()
    
#     # configure plot
#     plt.style.use("fast")
#     fig, axes = plt.subplots(
#         nrows=len(row_names),
#         ncols=len(col_names),
#         figsize=(40,40),
#         subplot_kw=dict(projection="3d")
#     )
#     for nth_row in range(len(row_names)):
#         for nth_col in range(len(col_names)):
#             # plot hue = sample
#             ax = axes[nth_row, nth_col]
#             r_name = row_names[nth_row]
#             c_name = col_names[nth_col]
#             subset = merged_data[(merged_data[row] == r_name) & (merged_data[col] == c_name)]
#             ax.view_init(elev=elev, azim=azim, roll=roll)
#             if subset.shape[0] == 0:
#                 continue
#             ax.scatter(
#                 subset["UMAP-1"],
#                 subset["UMAP-2"],
#                 subset["UMAP-3"],
#                 color=subset["Cluster"].map(palette),
#                 s=1
#             )
#             ax.set_xlabel("UMAP-1")
#             ax.set_ylabel("UMAP-2")
#             ax.set_zlabel("UMAP-3")
#             ax.set_title("{row}={r_name} | {col}={c_name} | n_cells={n} | {elev},{azim},{roll}".format(
#                 row=row,
#                 r_name=r_name,
#                 col=col,
#                 c_name=c_name,
#                 n=subset.shape[0],
#                 elev=elev,
#                 azim=azim,
#                 roll=roll
#             ))
#             ax.xaxis.pane.fill = False
#             ax.yaxis.pane.fill = False
#             ax.zaxis.pane.fill = False

#             # set color to white (or whatever is "invisible")
#             ax.xaxis.pane.set_edgecolor('w')
#             ax.yaxis.pane.set_edgecolor('w')
#             ax.zaxis.pane.set_edgecolor('w')

#             # get rid of the grid as well:
#             ax.grid(dashes=(1,1))
            
#     # save plot
#     plt.savefig(
#         output_prefix.with_suffix("").__str__() + "_hue=cluster.png",
#         dpi=300
#     )
    