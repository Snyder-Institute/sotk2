#!/usr/bin/env Rscript
# Then, you can merge the NMF outputs across ranks using the mergeNMFObjs function in sotk2.

# Do not run

library(NMF)
library(sotk2)

input_dir <- "/path/to/input_dir"
output_dir <- "/path/to/output_dir"

# Load NMF results and combine them
nmf_rank2 <- readRDS(file.path(input_dir, "NMF_rank_2.RDS"))
nmf_rank3 <- readRDS(file.path(input_dir, "NMF_rank_3.RDS"))
nmf_rank4 <- readRDS(file.path(input_dir, "NMF_rank_4.RDS"))

nmf_res <- mergeNMFObjs(list(`2` = nmf_rank2, `3` = nmf_rank3, `4` = nmf_rank4))
saveRDS(nmf_res, file.path(output_dir, "NMF_outputs.RDS"))

q('no')
