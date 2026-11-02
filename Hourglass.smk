from utils import functions
# N_PERMUTATIONS = 9999
# N_BOOTSTRAPS = 1000
N_PERMUTATIONS = 99
N_BOOTSTRAPS = 100

def ref_query_map(sample):
    mapping = {
        "egr_tenx_batch1": "egr_nextseq_batch1",
        "egr_nextseq_batch1": "egr_nextseq_batch1",
        "egr_nextseq_batch2": "egr_nextseq_batch1",
        "lch_tenx_batch1": "lch_tenx_batch1",
        "ptr_tenx_batch1": "ptr_tenx_batch1",
        "ptr_tenx_batch2": "ptr_tenx_batch1",
        "ptr_tenx_tsv1": "ptr_tenx_batch1",
        "ptr_tenx_tsv2": "ptr_tenx_batch1",
        "ptr_tenx_tsv3": "ptr_tenx_batch1",
        "ptr_tenx_tsv4": "ptr_tenx_batch1",
        "ptr_tenx_tsv5": "ptr_tenx_batch1",
        "tar_nextseq_batch1": "tar_nextseq_batch1"
    }
    try:
        return mapping[sample]
    except KeyError:
        raise ValueError(f"Unknown sample: {sample}")

rule transfer_cell_label:
    container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
    threads: 1
    resources:
        mem_mb=20000,
        runtime=10
    input:
        rmd="notebooks/Transfer_cell_label.Rmd",
        cell_label_map=lambda wildcards: "config/cell_label/%s.csv" % (ref_query_map(wildcards.sample)),
        path_10x_reference=lambda wildcards: functions.query(config["samples"], "name", ref_query_map(wildcards.sample))["results_old"]["10x"],
        path_10x_query=lambda wildcards: functions.query(config["samples"], "name", wildcards.sample)["results_old"]["10x"],
    output:
        path_outdir=directory("outputs/Hourglass/transfer_cell_label/{sample}/"),
        path_html="outputs/Hourglass/transfer_cell_label/{sample}/transfer_cell_label_{sample}.html",
        path_cell_label="outputs/Hourglass/transfer_cell_label/{sample}/cell_label.csv"
    params:
        reference_sample_id=lambda wildcards: ref_query_map(wildcards.sample),
        query_sample_id="{sample}",
        reference_cluster_names=lambda wildcards: f"config/cell_label/{ref_query_map(wildcards.sample)}.csv",
        n_pcs=30,
        resolution=0.8,
        prediction_score_threshold=0.50,
        prediction_margin_threshold=0.10,
        label_mode="cell",
        path_intermediate="outputs/Hourglass/transfer_cell_label/{sample}/rmd_intermediate/"
    log:
        "logs/Hourglass/transfer_cell_label/{sample}.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            '{input.rmd}',
            output_dir = '{output.path_outdir}',
            output_file = basename('{output.path_html}'),
            intermediates_dir = '{params.path_intermediate}',
            knit_root_dir = '~',
            params = list(
                reference_sample_id = '{params.reference_sample_id}',
                query_sample_id = '{params.query_sample_id}',
                path_10x_reference = '{input.path_10x_reference}',
                path_10x_query = '{input.path_10x_query}',
                reference_cluster_names = '{params.reference_cluster_names}',
                n_pcs = {params.n_pcs},
                resolution = {params.resolution},
                prediction_score_threshold = {params.prediction_score_threshold},
                prediction_margin_threshold = {params.prediction_margin_threshold},
                label_mode = '{params.label_mode}',
                path_cell_label = '{output.path_cell_label}'
            )
        )" > {log} 2>&1
        """

# def cell_label_samples(wildcards):
#     try:
#         samples = {
#             "ptr": [
#                 "ptr_tenx_batch1", # reference
#                 "ptr_tenx_batch2",
#                 "ptr_tenx_tsv1",
#                 "ptr_tenx_tsv2",
#                 "ptr_tenx_tsv3",
#                 "ptr_tenx_tsv4",
#                 "ptr_tenx_tsv5"
#             ],
#             "egr": {
#                 "egr_tenx_batch1", # reference 
#                 "egr_nextseq_batch1",
#                 "egr_nextseq_batch2"
#             }
#         }[wildcards.species]
#         return {
#             "samples": samples,
#             "path_10x": [functions.query(config["samples"], "name", s)["results_old"]["10x"] for s in samples]
#         }
        
#     except KeyError:
#         raise ValueError(f"No corresponding samples were defined for species: {wildcards.species}")

# rule cell_labelling:
#     container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
#     threads: 8
#     resources:
#         mem_mb=40,
#         runtime=120
#     input:
#         rmd="notebooks/_Consensus_cell_labelling_{species}.Rmd",
#         path_10x=lambda wildcards: cell_label_samples(wildcards)["path_10x"],
#         path_reference_label_mapping=lambda wildcards: 
#             "outputs/Hourglass/transfer_cell_label/%s/cell_label.csv" % (
#                 cell_label_samples(wildcards)["samples"][0]
#             )
#     output:
#         sample_manifest=temp("outputs/Hourglass/cell_labelling/{species}/sample-manifest.csv")
#     params:
#         samples=lambda wildcards: cell_label_samples(wildcards)["samples"]
#     log:
#         "logs/Hourglass/cell_labelling/{species}.log"
#     shell:
#         """
#         # write sample-manifest
#         mkdir -p $(dirname {output.sample_manifest})
#         echo Writing sample manifest to {output.sample_manifest} > {log}

#         samples=({params.samples})
#         path_10x=({input.path_10x})
#         echo 'sample,path_10x' > {output.sample_manifest}
#         for i in "${{!samples[@]}}"; do
#             echo "${{samples[i]}},${{path_10x[i]}}" >> "{output.sample_manifest}" 2>> {log}
#         done
#         """

rule intra_species_hourglass:
    container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
    threads: 76
    resources:
        mem_mb=400000,
        runtime=1200
    input:
        rmd="notebooks/Intra-species-hourglass.Rmd",
        path_10x=lambda wildcards: functions.query(config["samples"], "name", wildcards.sample_branched)["results_old"]["10x"],
        path_cell_label="outputs/Hourglass/transfer_cell_label/{sample_branched}/cell_label.csv",
    output:
        outdir=directory("outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}"),
        html="outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}/{sample_branched}-{metric}.html",
        rdata="outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}/{sample_branched}-{metric}.rdata",
        key_data="outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}/{sample_branched}-{metric}-key_data.csv",
        key_summary="outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}/{sample_branched}-{metric}-key_summary.csv"
    params:
        individual_id="{sample_branched}", # currently, individual_id is the same as sample name.
        condition="control", # currently, condition is set to control for all samples.
        batch=lambda wildcards: wildcards.sample_branched.split("_")[-1], # currently, batch is set to the last part of the sample name.
        keep_only_transfer_pass="TRUE", # currently, only keep cells that pass the transfer cell label step.
        n_bins=10,
        n_permutation=N_PERMUTATIONS,
        n_bootstrap=N_BOOTSTRAPS,
        window_dtw=0.4,
        metric="{metric}",
        path_intermediate="outputs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}/rmd_intermediate/"
    log:
        "logs/Hourglass/intra_species_hourglass/{sample_branched}/{metric}.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            '{input.rmd}',
            output_dir = '{output.outdir}',
            output_file = basename('{output.html}'),
            intermediates_dir = '{params.path_intermediate}',
            knit_root_dir = '~',
            clean = TRUE,
            envir = new.env(parent = globalenv()),
            params = list(
                sample_id = '{wildcards.sample_branched}',
                individual_id = '{params.individual_id}',
                condition = '{params.condition}',
                batch = '{params.batch}',
                path_10x = '{input.path_10x}',
                path_cell_label = '{input.path_cell_label}',
                keep_only_transfer_pass = as.logical(toupper('{params.keep_only_transfer_pass}')),
                n_bins = {params.n_bins},
                n_permutation = {params.n_permutation},
                n_bootstrap = {params.n_bootstrap},
                n_parallel = {threads},
                window_dtw = {params.window_dtw},
                metric = '{params.metric}',
                path_key_data = '{output.key_data}',
                path_key_summary = '{output.key_summary}',
                path_image = '{output.rdata}'
            )
        )" > {log} 2>&1 
        """

rule intra_species_hourglass_qi:
    container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
    threads: 76
    resources:
        mem_mb=400000,
        runtime=1200
    input:
        rmd="notebooks/Intra-species-hourglass_Qi.Rmd",
        path_rdata_in="rawdata/Qi_spatial/stall_fin_cnameChecked.Rdata",
        path_ray_middle_markers="test/Hourglass/Consensus_cell_labelling_ptr/Ray Precursor.csv",
        path_orthogroup="outputs/ortholog_info.csv"
    output:
        outdir=directory("outputs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}"),
        html="outputs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}/{sample_qi}-{metric}.html",
        key_data="outputs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}/{sample_qi}-{metric}-key_data.csv",
        key_summary="outputs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}/{sample_qi}-{metric}-key_summary.csv"
    params:
        individual_id="{sample_qi}", # currently, individual_id is the same as sample name.
        condition="control", # currently, all samples were from control conditions.
        batch=lambda wildcards: wildcards.sample_qi.split("_")[-1], # currently, batch is set to the last part of the sample name.
        assay="RNA",
        n_bins=10,
        n_permutation=N_PERMUTATIONS,
        n_bootstrap=N_BOOTSTRAPS,
        window_dtw=0.2,
        metric="{metric}",
        path_intermediate="outputs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}/rmd_intermediate/"
    log:
        "logs/Hourglass/intra_species_hourglass/{sample_qi}/{metric}.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            '{input.rmd}',
            output_dir = '{output.outdir}',
            output_file = basename('{output.html}'),
            intermediates_dir = '{params.path_intermediate}',
            knit_root_dir = '~',
            clean = TRUE,
            envir = new.env(parent = globalenv()),
            params = list(
                sample_id = '{wildcards.sample_qi}',
                biorep_id = '{params.batch}',
                individual_id = '{params.individual_id}',
                condition = '{params.condition}',
                batch = '{params.batch}',
                assay = '{params.assay}',
                path_rdata_in = '{input.path_rdata_in}',
                path_ray_middle_markers = '{input.path_ray_middle_markers}',
                path_orthogroup = '{input.path_orthogroup}',
                path_outdir = '{output.outdir}',
                n_bins = {params.n_bins},
                n_permutation = {params.n_permutation},
                n_bootstrap = {params.n_bootstrap},
                n_parallel = {threads},
                window_dtw = {params.window_dtw},
                metric = '{params.metric}',
                path_key_data = '{output.key_data}',
                path_key_summary = '{output.key_summary}'
            )
        )" > {log} 2>&1 
        """

rule inter_species_hourglass:
    container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
    threads: 60
    resources:
        mem_mb=400000,
        runtime=1200
    input:
        rmd="notebooks/Inter-species-hourglass (1).Rmd",
        path_10x_1=lambda wildcards: functions.query(config["samples"], "name", wildcards.sample1)["results_old"]["10x"],
        path_10x_2=lambda wildcards: functions.query(config["samples"], "name", wildcards.sample2)["results_old"]["10x"],
        path_cell_label_1="outputs/Hourglass/transfer_cell_label/{sample1}/cell_label.csv",
        path_cell_label_2="outputs/Hourglass/transfer_cell_label/{sample2}/cell_label.csv",
        orthogroup_table="outputs/ortholog_info.csv"
    output:
        outdir=directory("outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}"),
        html="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/{sample1}-{sample2}-{metric}.html",
        rdata="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/{sample1}-{sample2}-{metric}.rdata",
        key_data="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/{sample1}-{sample2}-{metric}-key_data.csv",
        key_summary="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/{sample1}-{sample2}-{metric}-key_summary.csv",
        orthogroups="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/{sample1}-{sample2}-{metric}-key_summary-orthogroups.csv"
    params:
        individual_id_1="{sample1}",
        individual_id_2="{sample2}",
        condition_1="control",
        condition_2="control",
        batch_1=lambda wildcards: wildcards.sample1.split("_")[-1],
        batch_2=lambda wildcards: wildcards.sample2.split("_")[-1],
        keep_only_transfer_pass="TRUE",
        n_bins=10,
        n_permutation=N_PERMUTATIONS,
        n_bootstrap=N_BOOTSTRAPS,
        window_dtw=0.1,
        metric="{metric}",
        path_intermediate="outputs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}/rmd_intermediate/"
    log:
        "logs/Hourglass/inter_species_hourglass/{sample1}/{sample2}/{metric}.log"
    shell:
        """
        Rscript -e "
        project_dir <- getwd()
        rmarkdown::render(
            '{input.rmd}',
            output_dir = '{output.outdir}',
            output_file = basename('{output.html}'),
            intermediates_dir = '{params.path_intermediate}',
            knit_root_dir = project_dir,
            clean = TRUE,
            envir = new.env(parent = globalenv()),
            params = list(
                sample_id_1 = '{wildcards.sample1}',
                sample_id_2 = '{wildcards.sample2}',
                individual_id_1 = '{params.individual_id_1}',
                individual_id_2 = '{params.individual_id_2}',
                condition_1 = '{params.condition_1}',
                condition_2 = '{params.condition_2}',
                batch_1 = '{params.batch_1}',
                batch_2 = '{params.batch_2}',
                path_10x_1 = '{input.path_10x_1}',
                path_10x_2 = '{input.path_10x_2}',
                path_cell_label_1 = '{input.path_cell_label_1}',
                path_cell_label_2 = '{input.path_cell_label_2}',
                path_orthogroups = '{input.orthogroup_table}',
                keep_only_transfer_pass = as.logical(toupper('{params.keep_only_transfer_pass}')),
                n_bins = {params.n_bins},
                n_permutation = {params.n_permutation},
                n_bootstrap = {params.n_bootstrap},
                n_parallel = {threads},
                window_dtw = {params.window_dtw},
                metric = '{params.metric}',
                path_key_data = '{output.key_data}',
                path_key_summary = '{output.key_summary}',
                path_image = '{output.rdata}'
            )
        )" > {log:q} 2>&1
        """