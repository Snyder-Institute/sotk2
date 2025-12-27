#!/usr/bin/env Rscript
# Once you have downloaded the required demonstration files from Zenodo, you can run the script below.
# If you saved the files to a user-defined directory, uncomment and set the download_dir parameter.

# download_dir <- "/path/to/download"

if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        download_dir <- download_dir
} else {
        download_dir <- tools::R_user_dir("sotk2", "data")
}

library(sotk2)

## Load NMF/cNMF results/objects
glass   <- readRDS(file.path(download_dir, "nmfRes_GLASS.RDS"))
ivygap  <- readRDS(file.path(download_dir, "nmfRes_IVYGAP.RDS"))
heiland <- readRDS(file.path(download_dir, "nmfRes_HEILAND.RDS"))

# data import
dataL        <- list(GLASS = glass, IVYGAP = ivygap, HEILAND = heiland)
rankL        <- list(GLASS = c(3:10, 15, 20), IVYGAP = c(3:10, 15, 20), HEILAND = seq(5, 30, 5))
dataCol      <- c("GLASS" = "cyan3", "IVYGAP" = "chartreuse1", "HEILAND" = "magenta")
corMethod    <- "spearman"

# Concat W matrics and calculate pairwise correlations
soSet <- SOSet(
        NMFobjL = dataL, 
        NMFrankL = rankL, 
        dataCol = dataCol, 
        corMet = corMethod
)
soSet

# threshold values/settings
corrCoefThre <- 0.3
seed         <- 1234
niter        <- 1000
commWeight   <- 100
cohortWeight <- 10

# Generate a correlation network and apply community search
soObj <- SOTK(
        SOSet = soSet, 
        coefThre = corrCoefThre, 
        seed = seed, 
        niter = niter, 
        commWeight = commWeight, 
        cohortWeight = cohortWeight
)
soObj

# Save the object
saveRDS(soObj, file.path(download_dir, "soObj.RDS"))
message("soObj.RDS was created.")

q("no")
