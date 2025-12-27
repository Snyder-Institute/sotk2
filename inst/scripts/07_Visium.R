#!/usr/bin/env Rscript
# Users can also overlay the Visium slide image with sotk2 results. 
# The script below provides a visualization example for a single sample (UKF269_T) and 
# generates a two-panel figure: 
# (left) Seurat SNN (shared nearest neighbor) clustering and 
# (right) community annotations from sotk2, where each spot can be assigned to multiple communities.

# You'll need to download the Visium objects to proceed with the following script.
# If you already downloaded other files from Zenodo to a user-defined directory, 
# uncomment and set the download_dir parameter below.
# IMPORTANT: If download_dir is not set (or is set incorrectly), this script will download 
# the full demo dataset again to the default location, which may take a long time.

# download_dir <- "/path/to/download" # where the demo .RDS files are located

## I don't intall Seruat package any computers as I don't use the Seurat package. I will update the following script later on my old computer.

options(timeout = max(1800, getOption("timeout"))) # 30 minutes
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        paths <- download_demo_data(set = "full", download_dir = download_dir)
} else {
        message("ALERT::download_dir is not set (or invalid).\nThe full demo dataset will be downloaded to the default directory, which may take a long time.")
        paths <- download_demo_data(set = "full")
}

library(Seurat)
library(ggplot2)
library(gridExtra)

heiland <- readRDS(file.path(download_dir, "UKF269_T_spots.RDS"))
SO <- readRDS(paths[["UKF269_T_Visium.RDS"]])
SO <- Seurat::UpdateSeuratObject(SO)
SO <- Seurat::RenameCells(SO, add.cell.id = "UKF269_T_")

ravCol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(ravCol) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12") 

heiCol <- c("#A6761D", "#66A61E", "#1B9E77", "#E7298A", "magenta", "grey50", "#D95F02", "#7570B3", "#8DD3C7", "beige", "#ececec")
names(heiCol) <- c("9-10", "10", "5-10", "5", "2-5", "Exc", "2", "2-10", "2-9", "2", "0")


pdf(file.path(outDir, "15_Visium_slide_UKF269_T.pdf"), width=14, height=7)
panels <- vector("list", 2)

panels[[1]] <- SpatialDimPlot(SO, group.by = "seurat_clusters", label = TRUE, label.size = 3, pt.size.factor = 250) + theme(legend.position = "right") + labs(title = "UKF269_T, SNN") + scale_fill_manual(values = ravCol[levels(SO@meta.data$seurat_clusters)]) 

excluded <- c(); community <- c()
for (spotName in rownames(SO@meta.data)) {
        if (spotName %in% names(heiland)) {
                comm <- paste(sort(unlist(heiland[spotName])), collapse="-")
                if (str_count(comm, "-") > 1) comm <- 0
        } else {
                comm <- "Exc"
                excluded <- c(excluded, spotName)
        }
        community <- c(community, comm)
}

community <- factor(community)
SO@meta.data$community <- community

panels[[2]] <- SpatialDimPlot(SO, group.by = "community", label = TRUE, label.size = 3, pt.size.factor = 250) + theme(legend.position = "right") + labs(title = "UKF269_T, Community") + scale_fill_manual(values = heiCol[levels(community)]) 
G <- grid.arrange(grobs = panels, top = "HEILAND (Ravi et al.)", ncol = 2, nrow = 1)

print(G)
dev.off()

q("no")
