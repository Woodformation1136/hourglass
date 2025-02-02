import pandas as pd
import argparse
import json
import sys
import os

# init argument parser
parser = argparse.ArgumentParser(description="Get cell color")
parser.add_argument("--clusters_old", help="Old 'Barcode' to 'Cluster' csv")
parser.add_argument("--clusters_new", help="New 'Barcode' to 'Cluster' csv")
parser.add_argument("--colors_old", help="List of comma separated colors for each cluster")
parser.add_argument("--sample_name", help="Sample name prefix to add to barcode", required=False)
parser.add_argument("barcode2color_csv", help="Output new 'Barcode' to 'Color' csv")

# get command line arguments
args = parser.parse_args()
clusters_old = args.clusters_old
colors_old = args.colors_old.split(",")
clusters_new = args.clusters_new
barcode2color_csv = args.barcode2color_csv    


# load data
clusters_new = pd.read_csv(clusters_new)
clusters_old = pd.read_csv(clusters_old)
members_new = {
    cluster_n: set(
        clusters_new[clusters_new["Cluster"] == cluster_n]["Barcode"].values
    ) for cluster_n in clusters_new["Cluster"].unique()
}
members_old = {
    cluster_n: set(
        clusters_old[clusters_old["Cluster"] == cluster_n]["Barcode"].values
    ) for cluster_n in clusters_old["Cluster"].unique()

}

# calculate pairwise jaccard similarity
jaccard = {}
best_match = {}
for cluster_n_new, barcodes_new in members_new.items():
    last_jaccard = 0
    for cluster_n_old, barcodes_old in members_old.items():
        jaccard_index = len(barcodes_new & barcodes_old) / len(barcodes_new | barcodes_old)
        jaccard[(cluster_n_new, cluster_n_old)] = jaccard_index
        if jaccard_index > last_jaccard:
            last_jaccard = jaccard_index
            best_match[cluster_n_new] = cluster_n_old
    # if no match found, assign to the last cluster
    if last_jaccard == 0:
        best_match[cluster_n_new] = -1

jaccard = pd.Series(jaccard).unstack()
jaccard.columns = [f"old_{cluster_n}" for cluster_n in jaccard.columns]
jaccard.index = [f"new_{cluster_n}" for cluster_n in jaccard.index]
jaccard = jaccard.round(2)
# write jaccard similarity to stderr
print(jaccard.to_string(), file=sys.stderr)

# assign new clusters with old colors
colors_old = {k+1: v for k, v in enumerate(colors_old)}
colors_old[-1] = "#00FFFF"
colors_new = {
    cluster_n: colors_old[best_match[cluster_n]] for cluster_n in clusters_new["Cluster"].unique()
}


# get cell color
clusters_new["Color"] = clusters_new["Cluster"].map(colors_new)

# add sample name to barcode
if args.sample_name:
    clusters_new["Barcode"] = clusters_new["Barcode"].apply(lambda x: f"{args.sample_name}:{x}")

# save data
os.makedirs(os.path.dirname(barcode2color_csv), exist_ok=True)
clusters_new.to_csv(barcode2color_csv, index=False)
