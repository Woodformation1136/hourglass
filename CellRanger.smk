from typing import List
import random

configfile: "configs/config.json"

wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),
    random_seed="[0-9]+",
    gff_ext="gff|gff3"

# # master random seed
# master_random_seed = 42
# random.seed(master_random_seed)
# random_seeds = [random.randint(0, 2**32 - 1) for i in range(100)]

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

# get tools
cellranger=config["tools"]["cellranger"]
vartrix=config["tools"]["vartrix"]

# rules
SAMPLES = [i["name"] for i in config["samples"] if not i.get("#TODO", False)]
RANDOM_SEED=17
rule all_for_hourglass:
    input:
        expand(
            "outputs/cellranger/reanalyze/{sample}_rs{random_seed}/",
            sample=SAMPLES,
            random_seed=RANDOM_SEED
        )
        # expand(
        #     "outputs/pyplot_umap/{sample}_rs{random_seed}.png",
        #     sample=SAMPLES,
        #     random_seed=RANDOM_SEED
        # ),
        # expand(
        #     "outputs/distance_matrix/umap/{sample}_rs{random_seed}.csv",
        #     sample=SAMPLES,
        #     random_seed=RANDOM_SEED
        # ),
        # expand(
        #     "outputs/distance_matrix/pca/{sample}.csv",
        #     sample=SAMPLES,
        #     random_seed=RANDOM_SEED
        # )

# rule all_for_studying_umap_consistency:
#     input:
#         expand(
#             "outputs/pyplot_umap/{sample}_rs{random_seed}.png",
#             sample="tung_batch1",
#             random_seed=random_seeds + [17]
#         ),
#         expand(
#             "outputs/distance_matrix/umap/{sample}_rs{random_seed}.csv",
#             sample="tung_batch1",
#             random_seed=random_seeds + [17]
#         ),
#         expand(
#             "outputs/distance_matrix/pca/{sample}.csv",
#             sample="tung_batch1",
#         )

# rule pca_distance_matrix:
#     conda: "envs/bio-pyps.yaml"
#     threads: 1
#     input:
#         "outputs/cellranger/reanalyze/{sample}_rs17/"
#     output:
#         "outputs/distance_matrix/pca/{sample}.csv"
#     params:
#         path_pca="outs/analysis/pca/gene_expression_10_components/projection.csv"
#     log:
#         "logs/distance_matrix/pca/{sample}.log"
#     shell:
#         """
#         python scripts/calculate_distance_matrix.py \
#             {input}/{params.path_pca} \
#             {output} \
#         2> {log} \
#         1> {log}
#         """

# rule umap_distance_matrix:
#     conda: "envs/bio-pyps.yaml"
#     threads: 1
#     input:
#         "outputs/cellranger/reanalyze/{sample}_rs{random_seed}/"
#     output:
#         "outputs/distance_matrix/umap/{sample}_rs{random_seed}.csv"
#     params:
#         path_umap_projection="outs/analysis/umap/gene_expression_2_components/projection.csv"
#     log:
#         "logs/distance_matrix/umap/{sample}_rs{random_seed}.log"
#     shell:
#         """
#         python scripts/calculate_distance_matrix.py \
#             {input}/{params.path_umap_projection} \
#             {output} \
#         2> {log} \
#         1> {log}
#         """

# rule pyplot_umap:
#     conda: "envs/bio-pyps.yaml"
#     threads: 1
#     input:
#         indir="outputs/cellranger/reanalyze/{sample}_rs{random_seed}/",
#         palette="configs/color_new.json"
#     output:
#         "outputs/pyplot_umap/{sample}_rs{random_seed}.png"
#     params:
#         path_umap_projection="outs/analysis/umap/gene_expression_2_components/projection.csv",
#         path_clusters="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv",
#     log:
#         "logs/pyplot_umap/{sample}_rs{random_seed}.log"
#     shell:
#         """
#         python scripts/plot_umap.py \
#             {input.indir}/{params.path_umap_projection} \
#             {params.path_clusters} \
#             {input.palette} \
#             {output} \
#         2> {log} \
#         1> {log}
#         """

# rule compare_ptr1_old_new_cluster_assignment:
#     threads: 1
#     input:
#         cluster_old="/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/cellranger_count_TenX_Ptr/outs/analysis/clustering/kmeans_10_clusters/clusters.csv",
#         cluster_new="outputs/cellranger/reanalyze/tung_batch1_rs17/outs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv",
#         color_old="configs/color_old.json"
#     output:
#         color_new="configs/color_new.json"
#     log:
#         "logs/compare_ptr1_old_new_cluster_assignment.log"
#     shell:
#         """
#         python scripts/compare_cluster_assignment.py \
#             {input.cluster_old} \
#             {input.cluster_new} \
#             {input.color_old} \
#             {output.color_new} \
#         2> {log} \
#         1> {log}
#         """

def cellranger_reanalyze_input(wildcards):
    sample_dict = query(config["samples"], "name", wildcards.sample)
    if sample_dict.get("matrix_path", None) is not None:
        return sample_dict["matrix_path"]
    else:
        return {
            "nextseq": "outputs/mars2hd5/%s.h5" % wildcards.sample,
            "10x": "outputs/cellranger/count/%s/outs/filtered_feature_bc_matrix.h5" % wildcards.sample
        }[query(config["samples"], "name", wildcards.sample)["platform"]]
    raise ValueError("No input found for %s" % wildcards.sample)

rule cellranger_reanalyze:
    threads: 8
    resources:
        mem_gb=80
    input:
        cellranger_reanalyze_input
    output:
        outdir=directory("outputs/cellranger/reanalyze/{sample}_rs{random_seed}"),
        tmpdir=temp(directory("{sample}_reanalyze_rs{random_seed}_tmp")),
        params=temp("configs/cellranger_reanalyze_{sample}_rs{random_seed}.param")
    params:
        num_principal_comps=10,
        max_clusters=20,
        umap_n_neighbors=50,
        umap_min_dist=0.5,
        random_seed=lambda wildcards: wildcards.random_seed
    log:
        "logs/cellranger/reanalysis/{sample}_rs{random_seed}.log"
    shell:
        """
        # create temporary parameter file
        echo "num_principal_comps,{params.num_principal_comps}" > {output.params}
        echo "max_clusters,{params.max_clusters}" >> {output.params}
        echo "umap_n_neighbors,{params.umap_n_neighbors}" >> {output.params}
        echo "umap_min_dist,{params.umap_min_dist}" >> {output.params}
        echo "random_seed,{params.random_seed}" >> {output.params}
        
        # run cellranger reanalyzes
        {cellranger} reanalyze \
            --id={output.tmpdir} \
            --params={output.params} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            --matrix={input} \
        2> {log} \
        1> {log} \
        && mkdir -p {output.outdir} \
        && mv {output.tmpdir}/* {output.outdir}
        """

def set_sample_names_param(wildcards):
    sample = query(config["samples"], "name", wildcards.sample)
    if "sample_names" in sample.keys():
        return "--sample=%s" % sample["sample_names"]
    else:
        return ""

rule cellranger_count:
    threads: 8
    resources:
        mem_gb=100
    input:
        transcriptome=lambda wildcards: 
            "outputs/cellranger/mkref/%s" % (
                query(config["samples"], "name", wildcards.sample)["species"],
            ),
        fastq_path=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["fastq_path"]
    output:
        outdir=directory("outputs/cellranger/count/{sample}"),
        out_matrix="outputs/cellranger/count/{sample}/outs/filtered_feature_bc_matrix.h5",
        tmpdir=temp(directory("{sample}_count_tmp"))
    params:
        create_bam="false",
        sample_names=set_sample_names_param,
        fastq_path=lambda wildcards: 
            ",".join(query(config["samples"], "name", wildcards.sample)["fastq_path"])
    log: 
        "logs/cellranger/count/{sample}.log"
    shell:
        """
        {cellranger} count \
            --id={output.tmpdir} \
            --transcriptome={input.transcriptome} \
            --fastqs={params.fastq_path} \
            --create-bam={params.create_bam} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            {params.sample_names} \
        2> {log} \
        1> {log} \
        && mkdir -p {output.outdir} \
        && mv {output.tmpdir}/* {output.outdir}
        """ 

rule mars2hd5:
    conda: "envs/bio-pyps.yaml"
    threads: 1
    input:
        umi_dir=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["umi_dir"],
        cds_list=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["cds_list"]
    output:
        h5="outputs/mars2hd5/{sample}.h5"
    params:
        umi_criteria=100
    log:
        "logs/mars2hd5/{sample}.log"
    shell:
        """
        python scripts/mars2hd5.py \
            --umi_dir {input.umi_dir} \
            --cds {input.cds_list} \
            --output_hd5 {output.h5} \
            --umi_criteria {params.umi_criteria} \
        2> {log} \
        1> {log}
        """

rule cellranger_mkref:
    threads: 16
    resources:
        mem_gb=64
    input:
        assembly=lambda wildcards: 
            query(config["references"], "species", wildcards.species)["assembly"],
        # gtf=lambda wildcards: 
        #     query(config["references"], "species", wildcards.species)["annotation-gtf"],
        gtf="outputs/cellranger/mkgtf/{species}.gtf-filtered"
    output:
        directory("outputs/cellranger/mkref/{species}")
    log:
        "logs/cellranger/mkref/{species}.log"
    shell:
        """
        {cellranger} mkref \
            --genome=$(basename {output}) \
            --fasta={input.assembly} \
            --genes={input.gtf} \
            --nthreads={threads} \
            --memgb={resources.mem_gb} \
        2> {log} \
        1> {log} \
        && mv $(basename {output}) $(dirname {output})
        """

# filter non-coding genes
rule cellranger_mkgtf:
    input:
        lambda wildcards: 
            query(config["references"], "species", wildcards.species)["annotation-gtf"],
    output:
        "outputs/cellranger/mkgtf/{species}.gtf-filtered"
    log:
        "logs/cellranger/mkgtf/{species}.log"
    shell:
        """
        {cellranger} mkgtf \
            {input} {output} \
            --attribute=gene_biotype:protein_coding \
        2> {log} \
        1> {log}             
        """

# convert gff to gtf
rule gffread_gff32gtf:
    input:
        gff="{any}.gff3"
    output:
        gtf="{any}.gtf"
    log:
        "{any}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
                gffread {input.gff} -T -o {output.gtf} \
        2> {log} \
        1> {log}
        """

# convert gff to gtf
rule gffread_gff2gtf:
    input:
        gff="{any}.gff"
    output:
        gtf="{any}.gtf"
    log:
        "{any}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
                gffread {input.gff} -T -o {output.gtf} \
        2> {log} \
        1> {log}
        """

rule add_gid_to_gtf:
    input:
        gtf="references/Lch/Lich.1_0.gtf"
    output:
        gtf="references/Lch/Lich.1_0_addedGID.gtf"
    log:
        "logs/add_gid_to_gff/lch.log"
    run:
        with open(input.gtf, "r") as f, open(output.gtf, "w") as g:
            for line in f:
                if line.startswith("#"):
                    g.write(line)
                else:
                    # add gene_id {id} if missin gene_id
                    fields = line.split("\t")
                    if not "gene_id" in fields[-1]:
                        extended_attribute = fields[-1].replace(
                            "transcript_id", "gene_id"
                        )
                        fields[-1] = fields[-1].strip() + " " + extended_attribute
                        
                    g.write("\t".join(fields))
