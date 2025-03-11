from typing import List
import random

configfile: "configs/config.json"

wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),
    random_seed="[0-9]+",
    gff_ext="gff|gff3",
    umap_rs="[0-9]+"


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
# ORTHO_INFO="/home/f06b22037/SSD2/JW/project_SingleCell/20220821_orthologous_14plus1/primary_transcripts/OrthoFinder/Results_Aug21/Orthogroups/Orthogroups_4Plus1_from14plus1_protein2gene.txt"
REFERENCE_SAMPLE="ptr_tenx_batch1"
CELLRANGER_RS=17
master_random_seed = 42
random.seed(master_random_seed)
random_seeds = [random.randint(0, 2**32 - 1) for i in range(100)]
UMAP_RS = [17] + random_seeds[:10]

# rule all
rule four_trees_batch1:
    input:
        inter_species=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_inter_species_{lineage}.svg",
            group=["four_trees_batch1"],
            cellranger_rs=CELLRANGER_RS,
            n=2,
            umap_rs=17,
            lineage=["fusiform", "ray"]
        ),
        slingshot_x=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_x.svg",
            group=["four_trees_batch1"],
            cellranger_rs=CELLRANGER_RS,
            n=2,
            umap_rs=17,
            sample=config["seurat_integrated"]["four_trees_batch1"]["member"]
        ),
        umap2d_x=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap2d/plot/rs{umap_rs}/x/",
            group=[
                "four_trees_batch1"
            ],
            cellranger_rs=CELLRANGER_RS,
            umap_rs=17
        ),
        slingshot_all=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_all.svg",
            group=["four_trees_batch1"],
            cellranger_rs=CELLRANGER_RS,
            n=2,
            umap_rs=17,
            sample=config["seurat_integrated"]["four_trees_batch1"]["member"]
        ),
        umap2d_all=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap2d/plot/rs{umap_rs}/all/",
            group=[
                "four_trees_batch1"
            ],
            cellranger_rs=CELLRANGER_RS,
            umap_rs=17
        )

rule four_trees:
    input:
        slingshot=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}.pdf",
            group=["four_trees"],
            cellranger_rs=CELLRANGER_RS,
            n=2,
            umap_rs=17,
            sample=config["seurat_integrated"]["four_trees"]["member"]
        ),
        umap2d=expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap2d/plot/rs{umap_rs}/",
            group=[
                "four_trees"
            ],
            cellranger_rs=CELLRANGER_RS,
            umap_rs=17
        )

# rule test_metacell:
#     input:
#         expand(
#             "test-metacell/{group}/rs{cellranger_rs}/{sample}_{param}.pdf",
#             group=["four_trees"],
#             cellranger_rs=CELLRANGER_RS,
#             sample=["egr_nextseq_batch1"],
#             param=[1,2,3]
#         )
# rule metacell:
#     input:
#         expand(
#             "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}.pdf",
#             group=["four_trees"],
#             cellranger_rs=CELLRANGER_RS,
#             sample=["egr_nextseq_batch1"]
#         )
    
# config["seurat_integrated"]["four_trees"]["member"]
rule umap3d:
    input:
        expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap3d/plot/rs{umap_rs}_hue=cluster.png",
            group=["four_trees"],
            cellranger_rs=CELLRANGER_RS,
            umap_rs=UMAP_RS
        )
rule umap2d:
    input:
        expand(
            "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap2d/plot/rs{umap_rs}_hue=cluster.png",
            group=[
                "four_trees_batch1"
            ],
            cellranger_rs=CELLRANGER_RS,
            umap_rs=17
        )

# # Metacell analysis (only for ptr_tenx_batch1)
# # ==============================================================================
# def get_test_metacell_params(wildcards):
#     param_dict = config["test-metacell"][wildcards.param]
#     params = ""
#     for k, v in param_dict.items():
#         params += f"--{k} {v} "
#     return params

# rule run_test_metacell:
#     input:
#         matrix_path="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrices/",
#         cell_ident="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
#     output:
#         output_pdf="test-metacell/{group}/rs{cellranger_rs}/{sample}_{param}.pdf"
#     params:
#         scdb_tmpdir="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}_{param}_tmp",
#         palette="configs/palette.json",
#         params=get_test_metacell_params,
#         sample_name=lambda wildcards: "%s_rs%s.csv" % (wildcards.sample, wildcards.cellranger_rs)
#     log:
#         "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}_{param}.log"
#     shell:
#         """
#         docker run \
#             {docker_mount} \
#             -u $(id -u) \
#             --rm \
#             chiaenu/metacell:0.3.41 \
#                 Rscript scripts/run-metacell/metacell.R \
#                     --matrix_path {input.matrix_path}/{params.sample_name} \
#                     --cell_ident {input.cell_ident} \
#                     --scdb_tmpdir {params.scdb_tmpdir} \
#                     --output_pdf {output.output_pdf} \
#                     --palette {params.palette} \
#                     {params.params} \
#             2> {log} 1> {log}
#         """

# def get_metacell_params(wildcards):
#     param_dict = config["seurat_integrated"][wildcards.group]["metacell_params"][wildcards.sample]
#     params = ""
#     for k, v in param_dict.items():
#         params += f"--{k} {v} "
#     return params

# rule run_metacell:
#     input:
#         matrix_path="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrices/",
#         cell_ident="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
#     output:
#         output_pdf="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}.pdf"
#     params:
#         scdb_tmpdir="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}_tmp",
#         palette="configs/palette.json",
#         params=get_metacell_params,
#         sample_name=lambda wildcards: "%s_rs%s.csv" % (wildcards.sample, wildcards.cellranger_rs)        
#     log:
#         "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/metacell/{sample}.log"
#     shell:
#         """
#         docker run \
#             {docker_mount} \
#             -u $(id -u) \
#             --rm \
#             chiaenu/metacell:0.3.41 \
#                 Rscript scripts/run-metacell/metacell.R \
#                     --matrix_path {input.matrix_path}/{params.sample_name} \
#                     --cell_ident {input.cell_ident} \
#                     --scdb_tmpdir {params.scdb_tmpdir} \
#                     --output_pdf {output.output_pdf} \
#                     --palette {params.palette} \
#                     {params.params} \
#             2> {log} 1> {log}
#         """


# Plot UMAP
# ==============================================================================
# "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
# "outputs/seurat_integration/groups/four_trees/rs{cellranger_rs}/cell_identities.csv"
# rule plot_umap_original:
#     conda: "envs/bio-pyps.yaml"
#     input:
#         umap_projection="outputs/cellranger/reanalyze/{sample}_rs{cellranger_rs}/outs/analysis/umap/gene_expression_2_components/projection.csv",
#         cell_identities="outputs/seurat_integration/groups/four_trees/rs{cellranger_rs}/cell_identities.csv"
#     output:
#         tmp_projection=temp(
#             "outputs/seurat_annotated/{sample}/cellranger{cellranger_rs}/umap{umap_rs}/proj.tmp"
#         ),
#         fpath="outputs/seurat_annotated/cellranger{cellranger_rs}/umap{umap_rs}/{sample}.png"
#     params:
#         palette="configs/palette.json"
#     log:
#         "logs/seurat_annotated/{sample}/cellranger{cellranger_rs}/umap{umap_rs}/plot_umap.log"
#     shell:
#         """
#         # prepending sample name to the first column of the projection file
#         bash scripts/prepend-sample-name.sh {input.umap_projection} {wildcards.sample}_rs{wildcards.cellranger_rs} {output.tmp_projection}

#         # plot umap
#         python scripts/plot-umap.py \
#             --cell_embedding {output.tmp_projection} \
#             --cell_identity {input.cell_identities} \
#             --output_fpath {output.fpath} \
#             --palette {params.palette} \
#         2> {log} 1> {log}
#         """

rule plot_umap_integrated_by_sample:
    conda:
        "envs/bio-pyps.yaml"
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv"
    output:
        umap_plot1="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}_hue=sample.png"
    params:
        palette="configs/palette.json",
        prefix="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}"
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/{group}.log"
    shell:
        """
        # s = 1
        python scripts/plot-umap-integrated-by-sample.py \
            --cell_embedding {input.umap_projection} \
            --output_prefix {params.prefix} \
            --palette {params.palette} \
            --n_components {wildcards.n} \
        2> {log} 1> {log}
        """

rule plot_umap_integrated_by_cluster_all:
    conda:
        "envs/bio-pyps.yaml"
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/cell_identities_unique_rs{umap_rs}.csv"
    output:
        umap_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}/all/")
    params:
        palette="configs/palette_all.json",
        prefix="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}/all/"
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/all/{group}.log"
    shell:
        """
        # RE-RUN
        python scripts/plot-umap-integrated.py \
            --cell_embedding {input.umap_projection} \
            --cell_identity {input.cell_identities} \
            --output_prefix {params.prefix} \
            --palette {params.palette} \
            --n_components {wildcards.n} \
        2> {log} 1> {log}
        """

rule plot_umap_integrated_by_cluster_x:
    conda:
        "envs/bio-pyps.yaml"
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/cell_identities_unique_rs{umap_rs}.csv"
    output:
        umap_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}/x/")
    params:
        palette="configs/palette.json",
        prefix="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/plot/rs{umap_rs}/x/"
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/x/{group}.log"
    shell:
        """
        # 
        python scripts/plot-umap-integrated.py \
            --cell_embedding {input.umap_projection} \
            --cell_identity {input.cell_identities} \
            --output_prefix {params.prefix} \
            --palette {params.palette} \
            --n_components {wildcards.n} \
        2> {log} 1> {log}
        """


# Run Slingshot
# ==============================================================================
rule plot_inter_species_ray:
    input:
        ptr="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_ptr_tenx_batch1_all/",
        egr="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_egr_nextseq_batch1_all/",
        tar="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_tar_nextseq_batch1_all/",
        lch="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_lch_tenx_batch1_all/",
        background_embedding="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv"
    output:
        output_plot="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_inter_species_ray.svg"
    params:
        palette="configs/palette_all.json",
        ratio=lambda _: 2.7 / 3.77
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/{group}/inter_species_ray.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/plot-multiple-lineages.R \
                    --lineages {input.ptr}/lineage_1.rds,{input.egr}/lineage_1.rds,{input.tar}/lineage_1.rds,{input.lch}/lineage_1.rds \
                    --background_embedding {input.background_embedding} \
                    --output_plot {output.output_plot} \
                    --ratio {params.ratio} \
            2> {log} 1> {log}
        """

rule plot_inter_species_fusiform:
    input:
        ptr="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_ptr_tenx_batch1_all/",
        egr="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_egr_nextseq_batch1_all/",
        lch="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_lch_tenx_batch1_all/",
        background_embedding="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv"
    output:
        output_plot="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_inter_species_fusiform.svg"
    params:
        palette="configs/palette_all.json",
        ratio=lambda _: 2.7 / 3.77
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/{group}/inter_species_fusiform.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/plot-multiple-lineages.R \
                    --lineages {input.ptr}/lineage_2.rds,{input.egr}/lineage_2.rds,{input.lch}/lineage_2.rds \
                    --background_embedding {input.background_embedding} \
                    --output_plot {output.output_plot} \
                    --ratio {params.ratio} \
            2> {log} 1> {log}
        """
        

        
# slingshot all
lineages_all = {
    "lch_tenx_batch1": "3,4,8/13,5,2,14/13,5,2,1",
    "egr_nextseq_batch1": "3,4,8/6,5,2,7/6,5,2,1",
    "ptr_tenx_batch1": "3,4,8/6,5,2,7/6,5,2,1",
    "tar_nextseq_batch1": "3,4,8/6,5,2,1"
}
rule run_slingshot_all:
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/cell_identities_unique_rs{umap_rs}.csv"
    output:
        sample_proj=temp("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}_{sample}_all.csv"),
        sample_cid=temp("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities_{n}d_{umap_rs}_{sample}_all.csv"),
        output_pdf="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_all.svg",
        output_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_all/")
    params:
        lineages=lambda wildcards: lineages_all[wildcards.sample],
        palette="configs/palette_all.json",
        sample=lambda wildcards: wildcards.sample,
        ratio=lambda _: 2.7 / 3.77
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/{group}/{sample}_all.log"
    shell:
        """
        # RE-RUN
        head -n1 {input.umap_projection} > {output.sample_proj}
        grep {params.sample} {input.umap_projection} >> {output.sample_proj}
        head -n1 {input.cell_identities} > {output.sample_cid}
        grep {params.sample} {input.cell_identities} >> {output.sample_cid}
        # re-run
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/slingshot.R \
                    --background_embedding {input.umap_projection} \
                    --umap_projection {output.sample_proj} \
                    --cell_identity {output.sample_cid} \
                    --lineages {params.lineages} \
                    --palette {params.palette} \
                    --output_pdf {output.output_pdf} \
                    --output_dir {output.output_dir} \
                    --ratio {params.ratio} \
            2> {log} 1> {log}
        """

# slingshot x part
lineages_x = {
    "lch_tenx_batch1": "3,4,8/13,5,2",
    "egr_nextseq_batch1": "3,4,8/6,5,2",
    "ptr_tenx_batch1": "3,4,8/6,5,2",
    "tar_nextseq_batch1": "3,4,8/6,5,2"
}
rule run_slingshot_x:
    input:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/cell_identities_unique_rs{umap_rs}.csv"
    output:
        sample_proj=temp("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}_{sample}_x.csv"),
        sample_cid=temp("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities_{n}d_{umap_rs}_{sample}_x.csv"),
        output_pdf="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_x.svg",
        output_curve1="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_x/lineage_1.rds",
        output_curve2="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_x/lineage_2.rds",
        output_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/slingshot/rs{umap_rs}_{sample}_x/")
    params:
        lineages=lambda wildcards: lineages_x[wildcards.sample],
        palette="configs/palette.json",
        sample=lambda wildcards: wildcards.sample,
        ratio=lambda _: 2.7 / 3.77
    log:
        "logs/seurat_annotated/cellranger{cellranger_rs}/umap{n}d{umap_rs}/{group}/{sample}_x.log"
    shell:
        """
        # re-re-run
        head -n1 {input.umap_projection} > {output.sample_proj}
        grep {params.sample} {input.umap_projection} >> {output.sample_proj}
        head -n1 {input.cell_identities} > {output.sample_cid}
        grep {params.sample} {input.cell_identities} >> {output.sample_cid}
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/slingshot:2.4.0 \
                Rscript scripts/run-slingshot/slingshot.R \
                    --background_embedding {input.umap_projection} \
                    --umap_projection {output.sample_proj} \
                    --cell_identity {output.sample_cid} \
                    --lineages {params.lineages} \
                    --palette {params.palette} \
                    --output_pdf {output.output_pdf} \
                    --output_dir {output.output_dir} \
                    --ratio {params.ratio} \
            2> {log} 1> {log}
        """


# Get cell identity by nearest neighbor
# ==============================================================================
# filtered cells will be assign cluster 13 (fusiform organizer) and 14 (fiber cells)
# cluster 11 cells will be ignored in the following analysis
rule umap_unique_cells:
    conda: "envs/bio-pyps.yaml"
    input:
        umap_embedding="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv",
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv",
        reference_cell_identity="outputs/seurat_integration/barcode2color/%s_rs{cellranger_rs}.csv" % (REFERENCE_SAMPLE)
    output:
        ref_cells=temp("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/ref_cells_{umap_rs}.txt"),
        cell_identities_unique="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/cell_identities_unique_rs{umap_rs}.csv"
    params:
        thresh=.90
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/rs{umap_rs}/umap_unique_cells.log"
    shell:
        """
        # rerun
        cat {input.reference_cell_identity} \
            | cut -d, -f1 \
            | tail -n+2 \
            > {output.ref_cells}
        python scripts/umap-unique-cell-in-lch.py \
            --umap_embedding {input.umap_embedding} \
            --ref_cells {output.ref_cells} \
            --input_cell_ident {input.cell_identities} \
            --output_cell_ident {output.cell_identities_unique} \
            --thresh {params.thresh} \
        2> {log} 1> {log}
        """

# filtered cells will be assign cluster 12
# cluster 11 cells will be ignored in the following analysis
rule get_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        embeddings="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv",
        reference_cell_ident=lambda wildcards:
            "outputs/seurat_integration/groups/%s/barcode2color/%s_rs%s_curated.csv" % (
                wildcards.group, REFERENCE_SAMPLE, wildcards.cellranger_rs
            )
    output:
        cell_identities="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/cell_identities.csv"
    params:
        thresh=.90
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/get_cell_identity.log"
    shell:
        """
        # re-run
        python scripts/get-cell-ident-by-NN.py \
            --embeddings {input.embeddings} \
            --reference_cell_ident {input.reference_cell_ident} \
            --cell_ident_csv {output.cell_identities} \
            --thresh {params.thresh} \
        2> {log} 1> {log}
        """

# filtered cells will be assign cluster 11
rule curate_ref_cell_identity:
    conda: "envs/bio-pyps.yaml"
    input:
        pca_embedding="outputs/seurat_integration/groups/{group}/rs17/pca_projection.csv",
        umap_embedding="outputs/seurat_integration/groups/{group}/rs17/umap2d/projections/rs17.csv",
        reference_cell_identity="outputs/seurat_integration/barcode2color/{sample}_rs{cellranger_rs}.csv"
    output:
        output="outputs/seurat_integration/groups/{group}/barcode2color/{sample}_rs{cellranger_rs}_curated.csv"
    params:
        thresh=.80,
    log:
        "logs/seurat_integration/curate_ref_cell_identity/groups/{group}/{sample}_rs{cellranger_rs}.log"
    shell:
        """
        # this line is to trigger re-run
        python scripts/mask-cluster-outlier.py \
            --embedding1 {input.pca_embedding} \
            --embedding2 {input.umap_embedding} \
            --reference_cell_ident {input.reference_cell_identity} \
            --output_csv {output.output} \
            --thresh {params.thresh} \
        2> {log} 1> {log}
        """


# Get reference cell identity
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
        pca_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv",
    output:
        umap_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/umap{n}d/projections/rs{umap_rs}.csv"
    params:
        params="configs/umap.json",
        random_state=lambda wildcards: wildcards.umap_rs,
        n_component=lambda wildcards: wildcards.n
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/run_umap/umap{n}d/{umap_rs}.log"
    shell:
        """
        python scripts/run-umap.py \
            --projections_in {input.pca_projection} \
            --projections_out {output.umap_projection} \
            --params {params.params} \
            --random_state {params.random_state} \
            --n_component {params.n_component} \
        2> {log} 1> {log}
        """

# Run Seurat
# ==============================================================================
def getParams(wildcards):
    param=""
    for k, v in config["seurat_integrated"][wildcards.group]["params"].items():
        param += f"--{k} {v} "
    return param

# TODO: add resource requirements of 50 GB memory
rule seurat:
    threads: 16
    resources:
        mem_gb=50
    input:
        ortho_umi_matrices=lambda wildcards:        
            expand(
                "outputs/seurat_integration/ortho_umi_count/{sample}_rs{cellranger_rs}.csv",
                sample=config["seurat_integrated"][wildcards.group]["member"],
                cellranger_rs=wildcards.cellranger_rs
            )
    output:
        pca_projection="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/pca_projection.csv",
        seurat_combined_object="outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/seurat_combined_object.rds"
    params:
        params=getParams
    log:
        "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/seurat.log"
    shell:
        """
        # now use old Seurat pipeline
        docker run \
            {docker_mount} \
            --rm \
            -u $(id -u):$(id -g) \
            satijalab/seurat:5.0.0 \
                Rscript scripts/run-seurat/seurat_integrate.R \
                    --threads {threads} \
                    --samples2integrate $(echo {input.ortho_umi_matrices} | sed 's/ /,/g') \
                    {params.params} \
                    --output_pca_projection_csv {output.pca_projection} \
                    --output_seurat_combined_object_rds {output.seurat_combined_object} \
        2> {log} 1> {log}
        """

# # split Seurat corrected expression matrix by sample
# rule split_corrected_expression_matrix:
#     input:
#         "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrix.csv"
#     output:
#         out_dir=directory("outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrices/")
#         # out_csv=lambda wildcards: expand(
#         #     "outputs/seurat_integration/groups/{group}/rs{cellranger_rs}/corrected_expression_matrices/{sample}_rs{cellranger_rs}.csv",
#         #     sample=config["seurat_integrated"][wildcards.group]["member"]
#         # )
#     params:
#         samples=lambda wildcards: config["seurat_integrated"][wildcards.group]["member"],
#         random_seed=lambda wildcards: wildcards.cellranger_rs 
#     log:
#         "logs/seurat_integration/groups/{group}/rs{cellranger_rs}/split_corrected_expression_matrix.log"
#     shell:
#         """
#         mkdir -p {output.out_dir}
#         for sample in {params.samples}; do
#             # Add first column name to file without inserting new line
#             echo -n Barcode > {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
#             head -n 1 {input} >> {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
#             grep $sample {input}  >> {output.out_dir}/${{sample}}_rs{params.random_seed}.csv
#         done
#         """

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
def get_ortho_umi_count_params(wildcards):
    id_mapping = ""
    if query(config["samples"], "name", wildcards.sample).get("id_mapping", False):
        id_mapping = "--id_mapping %s" % query(config["samples"], "name", wildcards.sample).get("id_mapping")
    return id_mapping

rule get_ortho_umi_count:
    conda: "envs/bio-pyps.yaml"
    input:
        feature_bc_mat="outputs/cellranger/reanalyze/{sample}_rs{cellranger_rs}/outs/filtered_feature_bc_matrix/",
        ortho_info="outputs/ortholog_info.csv"
    output:
        bc_ortho_mat="outputs/seurat_integration/ortho_umi_count/{sample}_rs{cellranger_rs}.csv"
    params:
        id_mapping=get_ortho_umi_count_params
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

# Run OrthoFinder
# ==============================================================================
