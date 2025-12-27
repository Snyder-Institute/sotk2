#!/usr/bin/env Rscript
# Users can also overlay the Visium slide image with sotk2 results. 
# The script below provides a visualization example for a single sample (UKF269_T) and 
# generates a two-panel figure: 
# (left) Seurat SNN (shared nearest neighbor) clustering and 
# (right) community annotations from sotk2, where each spot can be assigned to multiple communities.

# You'll need to download the Visium objects to proceed with the following script.
# If you already downloaded other files from Zenodo to a user-defined directory, 
# uncomment and set the download_dir parameter below.

# download_dir <- "/path/to/download" # where the demo .RDS files are located
# IMPORTANT: If download_dir is not set (or is set incorrectly), this script will download 
# the full demo dataset again to the default location, which may take a long time.

.demo_data_manifest <- function(set = c("core", "full")) {
        set <- match.arg(set)

        manifest_full <- data.frame(
                file = c(
                        "annot_GLASS.RDS",      # core
                        "expr_GLASS.RDS",       # core
                        "nmfRes_GLASS.RDS",     # core
                        "nmfRes_HEILAND.RDS",   # core
                        "nmfRes_IVYGAP.RDS",    # core
                        "UKF269_T_spots.RDS",
                        "UKF269_T_Visium.RDS"
                ),
                url = c(
                        "https://zenodo.org/records/18063318/files/annot_GLASS.RDS",      #  25 KB
                        "https://zenodo.org/records/18063318/files/expr_GLASS.RDS",       #  34 MB
                        "https://zenodo.org/records/18063318/files/nmfRes_GLASS.RDS",     #  13 MB
                        "https://zenodo.org/records/18063318/files/nmfRes_HEILAND.RDS",   #  57 MB
                        "https://zenodo.org/records/18063318/files/nmfRes_IVYGAP.RDS",    #  14 MB
                        "https://zenodo.org/records/18063318/files/UKF269_T_spots.RDS",   #  20 KB
                        "https://zenodo.org/records/18063318/files/UKF269_T_Visium.RDS"   # 542 MB (due to Visium v1 slide image)
                ),
                md5 = c(
                        "484c9d3637912c683a685182d6303d15",
                        "f270145499efdc1da08aac39ae4f4c78",
                        "269dc92869c8efad5ae539fe10cee9cb",
                        "016612979539a409fdf74c0692d97ffa",
                        "a12613ef9a08fc3328169ed7cbe83229",
                        "d2312d01b19695f5ec58bdeef231ea24",
                        "fad3a540e1280bbd9ee68fffc8df78d6"
                ),
                stringsAsFactors = FALSE
        )

        if (set == "full") {
                manifest <- manifest_full
        } else {
                manifest <- manifest_full[!grepl("^UKF269_T_", manifest_full$file), , drop = FALSE]
        }

        return(manifest)
}

.verify_md5 <- function(path, expected_md5) {
        if (is.na(expected_md5) || expected_md5 == "") {
                return(TRUE)
        }
        actual <- tools::md5sum(path)
        isTRUE(unname(actual) == expected_md5)
}

download_demo_data <- function(
        set = c("core", "full"),
        download_dir = tools::R_user_dir("sotk2", "data"),
        overwrite = FALSE) {

        set <- match.arg(set)

        if (!dir.exists(download_dir)) {
                dir.create(download_dir, recursive = TRUE)
        }

        manifest <- .demo_data_manifest(set = set)
        paths <- file.path(download_dir, manifest$file)

        message("Checking demo data files... (set = ", set, ")")

        for (i in seq_len(nrow(manifest))) {
                needs_download <- overwrite || !file.exists(paths[i])

                if (!needs_download) {
                        if (!.verify_md5(paths[i], manifest$md5[i])) {
                                message(sprintf(
                                        "Checksum mismatch for %s; re-downloading.",
                                        manifest$file[i]
                                ))
                                needs_download <- TRUE
                        }
                }

                if (needs_download) {
                        message(sprintf("Downloading %s", manifest$file[i]))
                        download.file(
                                url = manifest$url[i],
                                destfile = paths[i],
                                mode = "wb",
                                quiet = TRUE
                        )

                        if (!.verify_md5(paths[i], manifest$md5[i])) {
                                stop(sprintf(
                                        "Checksum verification failed for %s",
                                        manifest$file[i]
                                ))
                        }
                }
        }

        stats::setNames(paths, manifest$file)
}

options(timeout = max(1800, getOption("timeout"))) # 30 minutes
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        paths <- download_demo_data(set = "full", download_dir = download_dir)
} else {
        message("ALERT::download_dir is not set (or invalid).\nThe full demo dataset will be downloaded to the default directory, which may take a long time.")
        paths <- download_demo_data(set = "full")
}

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

pdf(file.path(download_dir, "17_Visium_slide_UKF269_T.pdf"), width=14, height=7)
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
