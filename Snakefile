# 
configfile: "config/config.json"

# 
include: "Hourglass.smk"

# Load custom functions used by all workflows
from utils import functions

# set wildcard constraints
wildcard_constraints:
    sample="(lch|ptr|egr|tar)_(tenx|nextseq)_(batch|tsv)[0-9]",
    sample1="lch_(tenx|nextseq)_(batch|tsv)[0-9]",
    sample2="(ptr|egr|tar)_(tenx|nextseq)_(batch|tsv)[0-9]",
    sample_qi="ptr_bmkgene_(1-7|4-7|3-7)",
    intra_workflow="Intra-species-hourglass(_protoplasting-excluded|_protoplasting-regressed)?",
    inter_workflow="Inter-species-hourglass(_protoplasting-excluded|_protoplasting-regressed)?"

rule summary:
    input:
        expand(
            "outputs/Hourglass-intra/{workflow}/summary.html",
            workflow=[
                "Intra-species-hourglass",
                "Intra-species-hourglass_protoplasting-excluded",
            ]
        ),
        expand(
            "outputs/Hourglass-inter/{workflow}/summary.html",
            workflow=[
                "Inter-species-hourglass",
                "Inter-species-hourglass_protoplasting-excluded"
            ]
        ),
        "outputs/Hourglass-Qi/summary.html"

rule all:
    input:
        # intra-species comparisons
        expand(
            "outputs/Hourglass/{intra_workflow}/{sample}/{metric}", 
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
                # "tar_nextseq_batch1",
                # "tar_nextseq_batch2",
                # "tar_tenx_batch1"
            ], 
            metric=["spearman", "pearson", "Jensen-Shannon"],
            intra_workflow=[
                "Intra-species-hourglass",
                # "Intra-species-hourglass_protoplasting-excluded"
                # "Intra-species-hourglass_protoplasting-regressed"
            ]
        ),
        expand(
            "outputs/Hourglass/Intra-species-hourglass/{sample}/{metric}", 
            sample=[
                "ptr_bmkgene_1-7",
                "ptr_bmkgene_3-7",
                "ptr_bmkgene_4-7",
            ], 
            metric=["spearman", "pearson", "Jensen-Shannon"]
        ),
        # inter-species comparisons
        expand(
            "outputs/Hourglass/{inter_workflow}/{sample1}/{sample2}/{metric}", 
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
            metric=["spearman", "pearson", "Jensen-Shannon"],
            inter_workflow=[
                "Inter-species-hourglass",
                # "Inter-species-hourglass_protoplasting-excluded"
                # "Inter-species-hourglass_protoplasting-regressed"
            ]
        )
