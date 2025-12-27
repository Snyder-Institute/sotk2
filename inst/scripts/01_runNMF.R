#!/usr/bin/env Rscript
# Usage: Rscript 01_runNMF.R GLASS $PROJECT_FOLDER/input $PROJECT_FOLDER/output 10
# You can parallelize NMF runs by executing each rank as an independent job, enabling efficient use of limited computational resources.

# Do not run

# User inputs
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
input_dir <- args[2]
output_dir <- args[3]
rank <- as.integer(args[4])

# Parameter settings
method <- "brunet" # default 
nr <- 1000
seed <- 123456
options <- "tv2p64" # # track + verbose level 2 + 64 threads
result_file_name <- paste0("NMF_rank_", rank)

# Loading
expr <- readRDS(file.path(input_dir, paste0(data_name, ".RDS")))
dim(expr)

# Run NMF
library(NMF)

if(nrow(expr) > 2) {
        NMFres <- nmf(expr, rank=rank, method=method, nrun=nr, seed=seed, .options=options)
        saveRDS(NMFres, file.path(output_dir, paste0(result_file_name, ".RDS")))

        exprRand <- randomize(expr)
        NMFresRand <- nmf(exprRand, rank=rank, method=method, nrun=nr, seed=seed, .options=options)
        saveRDS(NMFresRand, file.path(output_dir, paste0(result_file_name, "_rand.RDS")))
} else {
        message("Check your input.")
}

q("no")
