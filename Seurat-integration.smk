from typing import List
import random

configfile: "configs/config.json"

wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),
    random_seed="[0-9]+",
    gff_ext="gff|gff3"


# query function
def query(d:List[dict], k:str, v:str) -> dict:
    return [x for x in d if x[k] == v][0]

# create docker mount points
docker_mount = ""
for volume in config["volumes"]:
    docker_mount += "-v %s:%s:%s " % (
        volume["host"],
        volume["container"],
        volume["mode"]
    )
    if volume["is_workspace"]:
        docker_mount += "-w %s " % volume["container"]

# global variables
ORTHO_INFO="/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Ortholog_analysis/all_group_long_convertedID.csv"
REFERENCE_SAMPLE="ptr_tenx_batch1"
CELLRANGER_RS=17
master_random_seed = 42
random.seed(master_random_seed)
random_seeds = [random.randint(0, 2**32 - 1) for i in range(100)]
UMAP_RS=[17] + random_seeds[:10]

# rule all
SAMPLES = [i["name"] for i in config["samples"]]
RANDOM_SEED=17
GROUPS = [
    "ptr_tenx_batch1",
    "ptr_tenx",
    "ptr_all",
    "four_trees"
]
[
    "test_all"
]
rule all:
    input:
        # expand(
        #     "outputs/seurat_integration/groups/{group}/rs{CELLRANGER_RS}/pca_projection.csv",
        #     group=GROUPS
        # ),
        # expand(
        #     "outputs/seurat_integration/groups/{group}/rs{CELLRANGER_RS}/cell_identities.csv",
        #     group=GROUPS
        # ),
        expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap/plots/rs{umap_rs}/",
            group=GROUPS,
            cellranger_rs=CELLRANGER_RS,
            umap_rs=UMAP_RS
        )


# Plot UMAP
# ==============================================================================
rule plot_umap:
    conda: "envs/bio-pyps.yaml"
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
    output:
        plot_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap/plots/rs{umap_rs}/")
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/plot_umap/rs{umap_rs}.log"
    shell:
        """
        python scripts/plot-umap.py \
            {input.umap_projection} \
            {input.cell_identities} \
            {output.plot_dir} \
        2> {log} 1> {log}
        """


# Get cell identity by nearest neighbor
# ==============================================================================
rule get_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        embeddings="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv",
        reference_cell_ident=lambda wildcards:
            "outputs/seurat_integration/barcode2color/%s_rs%s.csv" % (
                REFERENCE_SAMPLE, wildcards.cellranger_rs
            )
    output:
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/get_cell_identity.log"
    shell:
        """
        python -i scripts/get-cell-ident-by-NN.py \
            --embeddings {input.embeddings} \
            --reference_cell_ident {input.reference_cell_ident} \
            --cell_ident_csv {output.cell_identities} \
        2> {log} 1> {log}
        """

rule get_ref_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        clusters_old=lambda wildcards: query(config["samples"], "name", wildcards.sample)["results_old"]["clustering"],
        clusters_new="outputs/cellranger/reanalyze/{sample}_rs{cellranger_rs}/outs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv"
    output:
        barcode2color_new="outputs/seurat_integration/barcode2color/{sample}_rs{cellranger_rs}.csv"
    params:
        colors_old=lambda wildcards: ",".join(query(config["samples"], "name", wildcards.sample)["results_old"]["colors"]).__repr__(),
        sample_name=lambda wildcards: "%s_rs%s" % (wildcards.sample, wildcards.cellranger_rs)
    
    log:
        "logs/seurat_integration/get_ref_cell_identity/{sample}_rs{cellranger_rs}.log"
    shell:
        """
        python scripts/get-cell-color.py \
            --clusters_old {input.clusters_old} \
            --clusters_new {input.clusters_new} \
            --colors_old {params.colors_old} \
            --sample_name {params.sample_name} \
            {output.barcode2color_new} \
        2> {log} 1> {log}
        """


# Run UMAP
# ==============================================================================
rule run_umap:
    conda: 
        "envs/umap-learn.yaml"
    input:
        pca_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv"
    output:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap/projections/rs{umap_rs}.csv"
    params:
        params="configs/umap.json",
        random_state=lambda wildcards: wildcards.umap_rs
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/run_umap/{umap_rs}.log"
    shell:
        """
        python scripts/run-umap.py \
            {input.pca_projection} \
            {output.umap_projection} \
            {params.params} \
            {params.random_state} \
        2> {log} 1> {log}
        """

# Run Seurat
# ==============================================================================
def getParams(wildcards):
    param=""
    for k, v in config["seurat_integrated"][wildcards.group]["params"].items():
        param += f"--{k}={v} "
    return param

rule seurat:
    input:
        ortho_umi_matrices=lambda wildcards:        
            expand(
                "outputs/seurat_integration/ortho_umi_count/{sample}_rs{cellranger_rs}.csv",
                sample=config["seurat_integrated"][wildcards.group]["member"],
                cellranger_rs=wildcards.cellranger_rs
            )
    output:
        pca_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv",
        corrected_expression_matrix="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrix.csv"
    params:
        params=getParams
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/seurat.log"
    shell:
        """
        docker run \
            {docker_mount} \
            --rm \
            -u $(id -u):$(id -g) \
            satijalab/seurat:5.0.0 \
                Rscript scripts/run-seurat/seurat_integrate.R \
                    --samples2integrate $(echo {input.ortho_umi_matrices} | sed 's/ /,/g') \
                    {params.params} \
                    --output_pca_projection_csv {output.pca_projection} \
                    --output_corrected_expression_matrix_csv {output.corrected_expression_matrix} \
        2> {log} 1> {log}
        """

# split Seurat corrected expression matrix by sample
rule split_corrected_expression_matrix:
    input:
        "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrix.csv"
    output:
        out_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrices/")
    params:
        samples=lambda wildcards: config["seurat_integrated"][wildcards.group]["member"],
        random_seed=lambda wildcards: wildcards.cellranger_rs 
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/split_corrected_expression_matrix.log"
    shell:
        """
        mkdir -p {output.out_dir}
        for sample in {params.samples}; do
            # Add first column name to file without inserting new line
            echo -n Barcode > {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
            head -n 1 {input} >> {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
            grep $sample {input} | \
                sed -re 's,([a-z0-9]+_)+,,g' \
                >> {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
        done
        """

# Convert UMI counts to Ortholog UMI counts
# ==============================================================================
# get ortholog information
rule get_ortholog_info:
    input:
        ORTHO_INFO
    output:
        "outputs/ortholog_info.csv"
    shell:
        """
        cp {input} {output}
        """

# get (barcode, ortholog UMI count) CSV from cellranger reanalyze output
rule get_ortho_umi_count:
    conda: "envs/bio-pyps.yaml"
    input:
        feature_bc_mat="outputs/cellranger/reanalyze/{sample}_rs{cellranger_rs}/outs/filtered_feature_bc_matrix/",
        ortho_info="outputs/ortholog_info.csv"
    output:
        bc_ortho_mat="outputs/seurat_integration/ortho_umi_count/{sample}_rs{cellranger_rs}.csv"
    params:
        id_mapping=lambda wildcards: "--id_mapping %s" if query(config["samples"], "name", wildcards.sample).get("id_mapping", False) else ""
    log:
        "logs/seurat_integration/ortho_umi_count/{sample}_rs{cellranger_rs}.log"
    shell:
        """
        python scripts/get-ortho-umi.py \
            --feature_bc_mat {input.feature_bc_mat} \
            --ortho_info {input.ortho_info} \
            --bc_ortho_mat {output.bc_ortho_mat} \
            {params.id_mapping} \
        2> {log} 1> {log}
        """
# ==============================================================================