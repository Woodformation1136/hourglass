import argparse
import pandas as pd

# create argument parser
argparser = argparse.ArgumentParser(
    description = "Get orthologous UMIs from a UMI table"
)
argparser.add_argument(
    "--feature_bc_mat",
    type=str,
    help="Filtered feature-barcode matrix (directory)"
)
argparser.add_argument(
    "--ortho_info",
    type=str,
    help="Orthologous information (.csv)"
)
argparser.add_argument(
    "--bc_ortho_mat",
    type=str,
    help="Output orthologous barcode matrix (.csv)"
)
argparser.add_argument(
    "--id_mapping",
    type=str,
    required=False,
    help="Mapping ids of feature-barcode matrix and ortho-info csv"
)

# parse arguments
args = argparser.parse_args()

# read orthologous information
ortho_info = pd.read_csv(args.ortho_info)
if args.id_mapping:
    if "Ath" in args.id_mapping:
        # get gene id to gene name mapping
        id_mapping = dict()
        with open(args.id_mapping) as f:
            for line in f:
                id_mapping[line.split("\t")[0]] = line.split("\t")[1].strip()
      
        # convert protein id to gene name
        ortho_info = ortho_info[ortho_info["species_name"] == "ArT"]
        gene_id = [gene.split(".")[0] for gene in ortho_info["gene_name"]]
        ortho_info["gene_name"] = [id_mapping.get(gene, None) for gene in gene_id]
    if "Osa" in args.id_mapping:
        # get transcript id to gene id mapping
        id_mapping = dict()
        with open(args.id_mapping) as f:
            for line in f:
                id_mapping[line.split("\t")[0]] = line.split("\t")[1].strip()
        
        # convert id in ortho_info to gene id
        ortho_info = ortho_info[ortho_info["species_name"] == "OrS"]
        ortho_info["gene_name"] = [id_mapping.get(gene, None) for gene in ortho_info["gene_name"]]
        
        
        
# TODO: read feature-barcode matrix
barcodes = pd.read_csv(
    f"{args.feature_bc_mat}/barcodes.tsv.gz",
    header=None,
    sep="\t"
)
features = pd.read_csv(
    f"{args.feature_bc_mat}/features.tsv.gz",
    header=None,
    sep="\t"
)
# read values
values = pd.read_csv(
    f"{args.feature_bc_mat}/matrix.mtx.gz",
    skiprows=3,
    header=None,
    sep=" ",
    names=["feature", "barcode", "count"]
)
# map feature-number to feature name in values
values["feature"] = values["feature"].map(
    lambda feature_number: features[1][feature_number - 1]
)
# map feature-name to ortho-feature
species_name = ortho_info[ortho_info["gene_name"] == values["feature"][0]]["species_name"].values[0]
species_ortho_info = ortho_info[ortho_info["species_name"] == species_name]
gene2cluster =  dict(zip(species_ortho_info["gene_name"], species_ortho_info["cluster_id"]))
values["ortho_feature"] = values["feature"].map(
    lambda gene_name: gene2cluster.get(gene_name, None)
)
values = values.dropna()

# group by ["barcode", "ortho_feature"] and sum counts
values = values.groupby(["barcode", "ortho_feature"]).sum().reset_index()
values["ortho_feature"] = values["ortho_feature"].astype(int)
# values["ortho_feature"] = "Cluster_" + values["ortho_feature"].astype(str)
# create barcode-feature matrix

matrix = values.pivot(
    index="barcode",
    columns="ortho_feature",
    values="count"
)
matrix.fillna(0, inplace=True)
matrix.index = barcodes[0]
matrix.columns = "Cluster_" + matrix.columns.astype(str)
matrix.index.name = "Barcode"

# write to file (.csv)
matrix.to_csv(args.bc_ortho_mat) 
