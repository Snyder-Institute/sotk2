#!/usr/bin/env Rscript
# Users can also overlay the Visium slide image with sotk2 results.
# The script below provides a visualization example for a single sample (UKF269_T) and
# generates a two-panel figure:
# (left) Seurat SNN (shared nearest neighbor) clustering and
# (right) community annotations from sotk2, where each spot can be assigned to multiple communities.

# This script requires the "full" demo bundle (~660 MB total, including the 542 MB Visium v1 slide image).
# You must explicitly set download_dir before running.

# download_dir <- "/path/to/download" # where the demo .RDS files are located

library(sotk2)

if (!exists("download_dir") || !is.character(download_dir) ||
        length(download_dir) != 1 || !nzchar(download_dir)) {
        stop(
                "ERROR::'download_dir' is not set.\n",
                "  Set download_dir to a writable directory before sourcing this script,\n",
                "  for example: download_dir <- \"~/Desktop/sotk2_demo\"\n",
                "  The Visium bundle is ~660 MB; this guard prevents an accidental download."
        )
}

# Ensure the Visium files are present locally; download the full bundle if needed.
options(timeout = max(1800, getOption("timeout"))) # 30 minutes
paths <- download_demo_data(set = "full", download_dir = download_dir)

library(Seurat)
library(stringr)
library(ggplot2)
library(gridExtra)

visSpots <- readRDS(file.path(download_dir, "UKF269_T_spots.RDS")) # Community
seuratObj <- readRDS(file.path(download_dir, "UKF269_T_Visium.RDS"))
seuratObj <- Seurat::UpdateSeuratObject(seuratObj)
seuratObj <- Seurat::RenameCells(seuratObj, add.cell.id = "269_T_")

snnCol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(snnCol) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

commCol <- c("#A6761D", "#66A61E", "#1B9E77", "#E7298A", "magenta", "grey50", "#D95F02", "#7570B3", "#8DD3C7", "beige", "#ececec")
names(commCol) <- c("9-10", "10", "5-10", "5", "2-5", "Exc", "2", "2-10", "2-9", "2", "0")

pdf(file.path(download_dir, "18_Visium_slide_UKF269_T.pdf"), width=14, height=7)
panels <- vector("list", 2)

panels[[1]] <- SpatialDimPlot(
                seuratObj,
                group.by = "seurat_clusters",
                label = TRUE, label.size = 3,
                pt.size.factor = 250) +
        theme(legend.position = "right") +
        labs(title = "UKF269_T, SNN") +
        scale_fill_manual(values = snnCol[levels(seuratObj@meta.data$seurat_clusters)])

excluded <- c(); community <- c()
for (spotName in rownames(seuratObj@meta.data)) {
        if (spotName %in% names(visSpots)) {
                comm <- paste(sort(unlist(visSpots[spotName])), collapse="-")
                if (str_count(comm, "-") > 1) comm <- 0
        } else {
                comm <- "Exc"
                excluded <- c(excluded, spotName)
        }
        community <- c(community, comm)
}

community <- factor(community)
seuratObj@meta.data$community <- community

panels[[2]] <- SpatialDimPlot(seuratObj,
                group.by = "community",
                label = TRUE, label.size = 3,
                pt.size.factor = 250) +
        theme(legend.position = "right") +
        labs(title = "UKF269_T, Community") +
        scale_fill_manual(values = commCol[levels(community)])

G <- grid.arrange(grobs = panels, top = "HEILAND (Ravi et al.)", ncol = 2, nrow = 1)

print(G)
dev.off()

q("no")
