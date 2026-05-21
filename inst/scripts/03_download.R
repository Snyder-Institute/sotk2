#!/usr/bin/env Rscript
# To run the scripts in the inst/scripts folder, you must first download the demonstration data from Zenodo.
# The download is handled by sotk2::download_demo_data() (~120 MB for the core bundle).
# If you want to use a user-defined output directory, uncomment and set the download_dir parameter.

# download_dir <- "/path/to/download"

library(sotk2)

options(timeout = max(1200, getOption("timeout"))) # 20 minutes
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        paths <- download_demo_data(set = "core", download_dir = download_dir)
} else {
        paths <- download_demo_data(set = "core")
}

q("no")
