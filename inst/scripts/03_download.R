#!/usr/bin/env Rscript
# To run the scripts in the inst/scripts folder, you must first download the demonstration data from Zenodo. 
# The script below downloads the files to a local directory; uncomment and set download_dir to your preferred destination.

# download_dir <- "/path/to/download"

# Internal manifest describing all demo files
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

options(timeout = max(1200, getOption("timeout"))) # 20 minutes
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        paths <- download_demo_data(set = "core", download_dir = download_dir)
} else {
        paths <- download_demo_data(set = "core")
}

q("no")
