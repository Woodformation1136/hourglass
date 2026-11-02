# 
configfile: "config/config.json"

# 
include: "Hourglass.smk"

# Load custom functions used by all workflows
from utils import functions

# set wildcard constraints
wildcard_constraints:
    sample_branched="(lch|ptr|egr)_(tenx|nextseq)_.+",
    sample1="lch_(tenx|nextseq)_.+",
    sample_qi="ptr_bmkgene_(1-7|4-7|3-7)"

rule all:
    input:
        "outputs/Hourglass/intra_species_hourglass/done.txt",
        "outputs/Hourglass/inter_species_hourglass/done.txt"

rule all_intra:
    input:
        expand(
            "outputs/Hourglass/intra_species_hourglass/{sample}/{metric}", 
            sample=[
                "ptr_tenx_batch1",
                "ptr_tenx_batch2",
                "ptr_tenx_tsv1",
                "ptr_tenx_tsv2",
                "ptr_tenx_tsv3",
                "ptr_tenx_tsv4",
                "ptr_tenx_tsv5",
                "egr_nextseq_batch1",
                "egr_nextseq_batch2",
                "egr_tenx_batch1",
                "lch_tenx_batch1",
                "ptr_bmkgene_1-7",
                "ptr_bmkgene_3-7",
                "ptr_bmkgene_4-7",
                # "tar_nextseq_batch1"
            ], 
            metric=["spearman", "pearson", "Jensen-Shannon"]
        )
    output:
        touch("outputs/Hourglass/intra_species_hourglass/done.txt")

rule all_inter:
    input:
        expand(
            "outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}", 
            sample1=[
                "lch_tenx_batch1"
            ], 
            sample2=[
                "ptr_tenx_batch1",
                "ptr_tenx_batch2",
                "ptr_tenx_tsv1",
                "ptr_tenx_tsv2",
                "ptr_tenx_tsv3",
                "ptr_tenx_tsv4",
                "ptr_tenx_tsv5",
                "egr_nextseq_batch1",
                "egr_nextseq_batch2",
                "egr_tenx_batch1",
            ],
            metric=["spearman", "pearson", "Jensen-Shannon"]
        )
    output:
        touch("outputs/Hourglass/inter_species_hourglass/done.txt")