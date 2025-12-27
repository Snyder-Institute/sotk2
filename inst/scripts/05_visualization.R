#!/usr/bin/env Rscript
# Using the Spatial Omics Set, you can generate diagnostic plots and network plots with the sotk2 package using the scripts below. 
# If you want to use a user-defined output directory, uncomment and set the download_dir parameter.

# download_dir <- "/path/to/download" # where soObj.RDS is located

if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        download_dir <- download_dir
} else {
        download_dir <- tools::R_user_dir("sotk2", "data")
}

if (!dir.exists(download_dir)) {
        dir.create(download_dir, recursive = TRUE)
        message(download_dir, " created.")
}

library(sotk2)

if (file.exists(file.path(download_dir, "soObj.RDS"))) {
        soObj <- readRDS(file.path(download_dir, "soObj.RDS"))
} else {
        stop("Run 04_runSotk2.R or download soObj.RDS from Zenodo.")
}

# visualization parameters
nodeSize        <- 10
nodeLabelSize   <- 2
edgeAlpha       <- 0.8

# visualization
plotCorrDensity(soObj@SOSet, filename = file.path(download_dir, "01_Stats_correlation_density.pdf")) # NULL to stdout
plotNetwork(soObj, label = FALSE, annot = "cohort", edgeAlpha = edgeAlpha, weighted = FALSE, filename = file.path(download_dir, "02_Network_Unweighted.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = TRUE, annot = "cohort", edgeAlpha = edgeAlpha, weighted = FALSE, filename = file.path(download_dir, "03_Network_Unweighted_lbl.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize) 
plotNetwork(soObj, label = FALSE, annot = "community", edgeAlpha = edgeAlpha, weighted = FALSE, filename = file.path(download_dir, "04_Network_Unweighted_Comm.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = TRUE, annot = "community", edgeAlpha = edgeAlpha, weighted = FALSE, filename = file.path(download_dir, "05_Network_Unweighted_Comm_lbl.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = FALSE, annot = "community", edgeAlpha = edgeAlpha, weighted = TRUE, filename = file.path(download_dir, "06_Network_Community.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = TRUE, annot = "community", edgeAlpha = edgeAlpha, weighted = TRUE, filename = file.path(download_dir, "07_Network_Community_lbl.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = FALSE, annot = "cohort", edgeAlpha = edgeAlpha, weighted = TRUE, filename = file.path(download_dir, "08_Network_Weighted.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)
plotNetwork(soObj, label = TRUE, annot = "cohort", edgeAlpha = edgeAlpha, weighted = TRUE, filename = file.path(download_dir, "09_Network_Weighted_lbl.pdf"), vertexSize = nodeSize, vertexLabelCex = nodeLabelSize)

# Community network property, number of GEPs in each community per dataset
statComm(soObj, filename = file.path(download_dir, "10_Number_of_GEPs_in_each_community.pdf"))

# Community-level network, layout
plotCommNetwork(soObj, vertexInfo = NULL, filename = file.path(download_dir, "11_community_network_layout.pdf"))

q("no")
