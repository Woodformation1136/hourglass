configfile: "configs/config.json"

# imports 
from typing import List

# get docker mount points
docker_mount = ""
for volume in config["volumes"]:
    docker_mount += "-v %s:%s:%s " % (
        volume["host"],
        volume["container"],
        volume["mode"]
    )
    if volume["is_workspace"]:
        docker_mount += "-w %s " % volume["container"]

# query function
def query(d:List[dict], k:str, v:str) -> dict:
    return [x for x in d if x[k] == v][0]

# master rule all
rule all:
    input:
        umap_x="outputs/Populus/umap/x/ptr_tenx_batch1.pdf",
        umap_all="outputs/Populus/umap/all/ptr_tenx_batch1.pdf",
        slingshot_umap_x="outputs/Populus/slingshot/x/ptr_tenx_batch1.pdf",
        slingshot_umap_all="outputs/Populus/slingshot/all/ptr_tenx_batch1.pdf",
        metacell_plot="outputs/Populus/metacell/ptr_tenx_batch1.pdf"

# transfer cell identity from Genome Biol. 2023 to new results
# ==============================================================================
rule get_ptr_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        clusters_old=lambda wildcards: query(config["samples"], "name", wildcards.sample)["results_old"]["clustering"],
        clusters_new="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv"
    output:
        barcode2color_new="outputs/Populus/barcode2color/{sample}_rs17.csv"
    params:
        colors_old=lambda wildcards: ",".join(query(config["samples"], "name", wildcards.sample)["results_old"]["colors"]).__repr__()  
    log:
        "logs/Populus/get_ref_cell_identity/{sample}_rs17.log"
    shell:
        """
        python scripts/get-cell-color.py \
            --clusters_old {input.clusters_old} \
            --clusters_new {input.clusters_new} \
            --colors_old {params.colors_old} \
            {output.barcode2color_new} \
        2> {log} 1> {log}
        """

rule curate_ptr_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        pca_embedding="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/pca/gene_expression_10_components/projection.csv",
        umap_embedding="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/umap/gene_expression_2_components/projection.csv",
        reference_cell_identity="outputs/Populus/barcode2color/{sample}_rs17.csv"
    output:
        "outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    params:
        thresh=.95
    log:
        "logs/Populus/curate_ptr_cell_identity/{sample}_rs17.log"
    shell:
        """
        python scripts/mask-cluster-outlier.py \
            --embedding1 {input.pca_embedding} \
            --embedding2 {input.umap_embedding} \
            --reference_cell_ident {input.reference_cell_identity} \
            --output_csv {output} \
            --thresh {params.thresh} \
        2> {log} 1> {log}
        """



# run metacell on CellRanger filtered UMI counts
# ==============================================================================
def get_metacell_params(wildcards):
    param_dict = query(config["samples"], "name", wildcards.sample)["metacell_params"]
    params = ""
    for k, v in param_dict.items():
        params += f"--{k}={v} "
    return params

rule metacell:
    input:
        matrix_path="outputs/cellranger/reanalyze/{sample}_rs17/outs/filtered_feature_bc_matrix",
        cell_ident="outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    output:
        output_pdf="outputs/Populus/metacell/{sample}.pdf"
    params:
        scdb_tmpdir="outputs/Populus/metacell/{sample}_scdb_tmpdir",
        palette="configs/palette.json",
        params=get_metacell_params
    log:
        "logs/Populus/metacell/{sample}.log"
    shell:
        """
        # metacell
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/metacell:0.3.41 \
                Rscript scripts/run-metacell/metacell.R \
                    --matrix_path {input.matrix_path} \
                    --cell_ident {input.cell_ident} \
                    --scdb_tmpdir {params.scdb_tmpdir} \
                    --output_pdf {output.output_pdf} \
                    --palette {params.palette} \
                    {params.params} \
            2> {log} 1> {log}
        """

# run slingshot on cellranger reanalyzed umap
# ==============================================================================
rule slingshot_x:
    input:
        umap_projection="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/umap/gene_expression_2_components/projection.csv",
        cell_identity="outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    output:
        output_pdf="outputs/Populus/slingshot/x/{sample}.pdf"
    params:
        lineages="3,4,8/6,5,2",
        palette="configs/palette.json"
    log:    
        "logs/Populus/slingshot/{sample}_x.log"
    shell:
        """
        # rerun slingshot
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/slingshot.R \
                    --background_embedding {input.umap_projection} \
                    --umap_projection {input.umap_projection} \
                    --cell_identity {input.cell_identity} \
                    --lineages {params.lineages} \
                    --palette {params.palette} \
                    --output_pdf {output.output_pdf} \
            2> {log} 1> {log}
        """

rule slingshot_all:
    input:
        umap_projection="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/umap/gene_expression_2_components/projection.csv",
        cell_identity="outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    output:
        output_pdf="outputs/Populus/slingshot/all/{sample}.pdf"
    params:
        lineages="3,4,8/6,5,2,7/6,5,2,1",
        palette="configs/palette_all.json"
    log:    
        "logs/Populus/slingshot/{sample}_all.log"
    shell:
        """
        # re-run
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/slingshot.R \
                    --background_embedding {input.umap_projection} \
                    --umap_projection {input.umap_projection} \
                    --cell_identity {input.cell_identity} \
                    --lineages {params.lineages} \
                    --palette {params.palette} \
                    --output_pdf {output.output_pdf} \
            2> {log} 1> {log}
        """

# plot umap
# ==============================================================================
rule plot_umap_x:
    conda:
        "envs/bio-pyps.yaml"
    input:
        umap="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/umap/gene_expression_2_components/projection.csv",
        cell_identity="outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    output:
        "outputs/Populus/umap/x/{sample}.pdf"
    params:
        palette="configs/palette.json",
        prefix="outputs/Populus/umap/x/{sample}"
    log:
        "logs/Populus/umap/x/{sample}.log"
    shell:
        """
        # rerun
        python scripts/plot-umap-original.py \
            --cell_embedding {input.umap} \
            --cell_identity {input.cell_identity} \
            --palette {params.palette} \
            --output_fpath {output} \
        2> {log} 1> {log}
        """

rule plot_umap:
    conda:
        "envs/bio-pyps.yaml"
    input:
        umap="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/umap/gene_expression_2_components/projection.csv",
        cell_identity="outputs/Populus/barcode2color/{sample}_rs17_curated.csv"
    output:
        "outputs/Populus/umap/all/{sample}.pdf"
    params:
        palette="configs/palette_all.json",
        prefix="outputs/Populus/umap/all/{sample}"
    log:
        "logs/Populus/umap/all/{sample}.log"
    shell:
        """
        # rerun
        python scripts/plot-umap-original.py \
            --cell_embedding {input.umap} \
            --cell_identity {input.cell_identity} \
            --palette {params.palette} \
            --output_fpath {output} \
        2> {log} 1> {log}
        """