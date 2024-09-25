from typing import List
import random

configfile: "configs/config.json"

wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),
    random_seed="[0-9]+"

# master random seed
master_random_seed = 42
random.seed(master_random_seed)
random_seeds = [random.randint(0, 2**32 - 1) for i in range(100)]

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
rule all:
    input:
        expand(
            "outputs/pyplot_umap/{sample}_rs{random_seed}.png",
            sample=[i["name"] for i in config["samples"]],
            random_seed=random_seeds + [17]
        ),
        expand(
            "outputs/distance_matrix/umap/{sample}_rs{random_seed}.csv",
            sample=[i["name"] for i in config["samples"]],
            random_seed=random_seeds + [17]
        ),
        expand(
            "outputs/distance_matrix/pca/{sample}.csv",
            sample=[i["name"] for i in config["samples"]],
        )

rule pca_distance_matrix:
    conda: "envs/bio-pyps.yaml"
    threads: 1
    input:
        "outputs/cellranger/reanalyze/{sample}_rs17/"
    output:
        "outputs/distance_matrix/pca/{sample}.csv"
    params:
        path_pca="outs/analysis/pca/gene_expression_10_components/projection.csv"
    log:
        "logs/distance_matrix/pca/{sample}.log"
    shell:
        """
        python scripts/calculate_distance_matrix.py \
            {input}/{params.path_pca} \
            {output} \
        2> {log} \
        1> {log}
        """

rule umap_distance_matrix:
    conda: "envs/bio-pyps.yaml"
    threads: 1
    input:
        "outputs/cellranger/reanalyze/{sample}_rs{random_seed}/"
    output:
        "outputs/distance_matrix/umap/{sample}_rs{random_seed}.csv"
    params:
        path_umap_projection="outs/analysis/umap/gene_expression_2_components/projection.csv"
    log:
        "logs/distance_matrix/umap/{sample}_rs{random_seed}.log"
    shell:
        """
        python scripts/calculate_distance_matrix.py \
            {input}/{params.path_umap_projection} \
            {output} \
        2> {log} \
        1> {log}
        """

rule pyplot_umap:
    conda: "envs/bio-pyps.yaml"
    threads: 1
    input:
        "outputs/cellranger/reanalyze/{sample}_rs{random_seed}/"
    output:
        "outputs/pyplot_umap/{sample}_rs{random_seed}.png"
    params:
        path_umap_projection="outs/analysis/umap/gene_expression_2_components/projection.csv",
        path_clusters="outputs/cellranger/reanalyze/{sample}_rs17/outs/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv",
    log:
        "logs/pyplot_umap/{sample}_rs{random_seed}.log"
    shell:
        """
        python scripts/plot_umap.py \
            {input}/{params.path_umap_projection} \
            {params.path_clusters} \
            {output} \
        2> {log} \
        1> {log}
        """

rule cellranger_reanalyze:
    threads: 24
    resources:
        mem_gb=80
    input:
        "outputs/cellranger/count/{sample}"
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
            --matrix={input}/outs/filtered_feature_bc_matrix.h5 \
        2> {log} \
        1> {log} \
        && mkdir -p {output.outdir} \
        && mv {output.tmpdir}/* {output.outdir}
        """

rule cellranger_count:
    threads: 32
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
        tmpdir=temp(directory("{sample}_count_tmp"))
    params:
        create_bam="true"
    log: 
        "logs/cellranger/count/{sample}.log"
    shell:
        """
        {cellranger} count \
            --id={output.tmpdir} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.fastq_path} \
            --create-bam={params.create_bam} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
        2> {log} \
        1> {log} \
        && mkdir -p {output.outdir} \
        && mv {output.tmpdir}/* {output.outdir}
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
rule gffread_gff2gtf:
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