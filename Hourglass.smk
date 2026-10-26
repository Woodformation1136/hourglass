from utils import functions

rule intra_species_hourglass:
    threads: 76
    resources:
        mem_mb=190000,
        runtime=1200
    container: "docker://chiaenu/rstudio-seurat:5.0.0-spacexr-slingshot2.4.0-dtw-sccore-dosnow"
    input:
        rmd = "notebooks/Intra-species_hourglass.Rmd",
        cellranger_out = lambda wildcards: functions.query(config["samples"], "name", wildcards.sample)["results_old"]["10x"],
        cell_label = "config/cell_label/{sample}.csv",
    output:
        outdir = directory("outputs/Hourglass/{sample}/{metric}"),
        html = "outputs/Hourglass/{sample}/{metric}/{sample}-{metric}.html",
        rdata = "outputs/Hourglass/{sample}/{metric}/{sample}-{metric}.rdata"
    log:
        "logs/Hourglass/{sample}/{metric}.log"
    params:
        metric = "{metric}",
        n_bins = 10,
        n_permutation = 9999,
        n_bootstrap = 1000,
        window_dtw = 0.5
    shell:
        """
        Rscript -e "rmarkdown::render(
            '{input.rmd}',
            output_dir = '{output.outdir}',
            output_file = '{output.html}',
            knit_root_dir = '~',
            params = list(
                path_10x = '{input.cellranger_out}',
                cluster_names = '{input.cell_label}',
                n_bins = {params.n_bins},
                n_permutation = {params.n_permutation},
                n_bootstrap = {params.n_bootstrap},
                n_parallel = {threads},
                metric = '{params.metric}',
                window_dtw = {params.window_dtw},
                path_image = '{output.rdata}'
            )
        )" > {log} 2>&1 
        """