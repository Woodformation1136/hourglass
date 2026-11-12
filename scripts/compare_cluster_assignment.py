import pandas as pd
import json
import sys

argv = sys.argv[1:]
if len(argv) < 4:
    print("Usage: compare_cluster_assignment.py <clusters.old> <clusters.new> <colors.old> <colors.new>")
    sys.exit(1)

clusters_old, clusters_new, colors_old, output_new = argv

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

jaccard = pd.Series(jaccard).unstack()
jaccard.columns = [f"old_{cluster_n}" for cluster_n in jaccard.columns]
jaccard.index = [f"new_{cluster_n}" for cluster_n in jaccard.index]
jaccard = jaccard.round(2)
print(jaccard, file=sys.stderr)

# assign new clusters with old colors
colors_old = json.load(open(colors_old))
colors_old = {int(k): v for k, v in colors_old.items()}
colors_new = {
    str(cluster_n): colors_old[best_match[cluster_n]] for cluster_n in clusters_new["Cluster"].unique()
}

# save new colors
with open(output_new, "w") as f:
    json.dump(colors_new, f)