library(slingshot)
library(rjson)
library(magrittr)
source("scripts/run-slingshot/utils.R")


# parse command line arguments and load data
parse_arguments()
background_embedding <- read.csv(background_embedding, row.names = 1) # Barcode, UMAP.1, UMAP.2
lineages <- strsplit(lineages, ",")[[1]] %>% lapply(readRDS) # list of list(curve_points, color)
ratio <- as.numeric(ratio)

# plot
width <- 7
height <- width * ratio
if (file.exists(output_plot)) {
    message("Removing existing plot file: ", output_plot)
    unlink(output_plot)
}
svg(output_plot, width = width, height = height)
par(mar = c(0, 0, 0, 0))
plot(
    background_embedding[, 1],
    background_embedding[, 2],
    col = "#D4D4D4",
    pch = 20,
    cex = 0.3,
    ann = FALSE,
    xaxt = "n",
    yaxt = "n",
    axes = FALSE
)
for (lineage in lineages) {
    points(
        lineage$curve_points,
        col=lineage$color,
        pch = 20,
        cex = 2
    )
}
dev.off()