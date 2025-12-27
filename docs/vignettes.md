# Introduction

  - This demo script illustrates the use of the `sotk2` package through a
  step-by-step workflow:
    1.  It demonstrates how to construct a Spatial Omics Set (`soSet`) by
      integrating (i) NMF-derived outputs as inputs and (ii) an
      all-pairs correlation matrix as the primary output.
    2.  The workflow then proceeds to create a spatial omics object
      (`soObj`), in which a metagene network is inferred by thresholding
      correlation strengths (for example, Spearman's rho \> 0.5).
    3.  A community detection algorithm is subsequently applied to the
      resulting network to identify metagene communities, which can be
      interpreted as candidate biological modules.
    4.  Following community detection, each community can be
      systematically annotated using sample-level metadata (for example,
      molecular subtype labels such as the Verhaak classification) to
      facilitate biological interpretation.
    5.  Users may also incorporate external annotations or
      project-specific features onto the network to support tailored
      interpretation and downstream analyses.

## Data download to run the demos

  - Demo data are available from Zenodo (<a href="https://doi.org/10.5281/zenodo.18063318" target="_blank" rel="noopener noreferrer">https://doi.org/10.5281/zenodo.18063318</a>).
  - You may download the files directly from the Zenodo record or use the
  provided helper script to retrieve them programmatically.
  - Two data bundles are provided: **core** and **full**
    1. The **core** bundle contains five files comprising cNMF outputs for
    the GLASS, IVYGAP, and HEILAND datasets, along with the GLASS cohort
    expression matrix and sample-level annotation metadata.
    2. The **full** bundle includes two additional files; however,
    `UKF269_T_Visium.RDS` (which contains Visium v1 slide image) is
    approximately 542 MB. Because of its size, we recommend downloading
    the core bundle first to begin running the demos, and downloading
    `UKF269_T_Visium.RDS` in the background. This Visium object is
    required for the visualization steps later in the workflow.

``` r
# The script below downloads the files to a local directory;
# uncomment and set download_dir to your preferred destination.

# download_dir <- "/path/to/download"

.demo_data_manifest <- function(set = c("core", "full")) {
        set <- match.arg(set)

        manifest_full <- data.frame(
                file = c(
                        "annot_GLASS.RDS",      # core,  25 KB
                        "expr_GLASS.RDS",       # core,  34 MB
                        "nmfRes_GLASS.RDS",     # core,  13 MB
                        "nmfRes_HEILAND.RDS",   # core,  57 MB
                        "nmfRes_IVYGAP.RDS",    # core,  14 MB
                        "UKF269_T_spots.RDS",   # full,  20 KB
                        "UKF269_T_Visium.RDS".  # full, 542 MB (Visium v1 slide image)
                ),
                url = c(
                        "https://zenodo.org/records/18063318/files/annot_GLASS.RDS",
                        "https://zenodo.org/records/18063318/files/expr_GLASS.RDS",
                        "https://zenodo.org/records/18063318/files/nmfRes_GLASS.RDS",
                        "https://zenodo.org/records/18063318/files/nmfRes_HEILAND.RDS",
                        "https://zenodo.org/records/18063318/files/nmfRes_IVYGAP.RDS",
                        "https://zenodo.org/records/18063318/files/UKF269_T_spots.RDS",
                        "https://zenodo.org/records/18063318/files/UKF269_T_Visium.RDS"
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
                manifest <- manifest_full[!grepl("^UKF269_T_", 
                        manifest_full$file), , drop = FALSE]
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
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 
        && nzchar(download_dir)) {
        paths <- download_demo_data(set = "core", download_dir = download_dir)
} else {
        paths <- download_demo_data(set = "core")
}
```

    ## Checking demo data files... (set = core)

## Directory settings

  - This block defines the working directory used throughout the demo to
  store downloaded inputs and generated outputs.
  - If a valid `download_dir` variable already exists in the current R
  session, the script reuses that user-specified location.
  - Otherwise, it defaults to a standard, package-specific data directory
  returned by `tools::R_user_dir("sotk2", "data")`, which provides a
  consistent and user-writable storage path across platforms.

``` r
if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 
        && nzchar(download_dir)) {
        download_dir <- download_dir
} else {
        download_dir <- tools::R_user_dir("sotk2", "data")
}
```

## Load the sotk2 library

  - If `sotk2` is not yet installed, please refer to the project
  repository for installation instructions: <a href="https://github.com/Snyder-Institute/sotk2" target="_blank" rel="noopener noreferrer">https://github.com/Snyder-Institute/sotk2</a>.
  - After loading `sotk2`, its required dependencies (including `igraph`,
  `methods`, `NMF`, `RColorBrewer`, and `stringr`) are loaded to support
  the full demonstration workflow.

``` r
# install.packages("devtools")
# devtools::install_github("Snyder-Institute/sotk2")

library(sotk2)
```

## Load cNMF objects

  - In this section, we load the precomputed cNMF result objects for three
  cohorts (GLASS, IVYGAP, and HEILAND) from the downloaded demo files.
  - These cohort-specific objects are then assembled into a single named
  list (`dataL`), which serves as the primary input for subsequent
  `sotk2` workflows.
  - We additionally define:
    1.  The candidate factorization ranks to consider for each cohort
      (`rankL`)
    2.  A cohort-specific color palette (`dataCol`) to ensure consistent
      visualization across figures
    3.  The correlation method (`corMethod`) used to quantify metagene
      similarity during the downstream integration and
      network-construction steps

``` r
glass     <- readRDS(file.path(download_dir, "nmfRes_GLASS.RDS"))
ivygap    <- readRDS(file.path(download_dir, "nmfRes_IVYGAP.RDS"))
heiland   <- readRDS(file.path(download_dir, "nmfRes_HEILAND.RDS"))

dataL     <- list(
                     GLASS = glass, 
                     IVYGAP = ivygap, 
                     HEILAND = heiland
             )
rankL     <- list(
                     GLASS = c(3:10, 15, 20), 
                     IVYGAP = c(3:10, 15, 20), 
                     HEILAND = seq(5, 30, 5)
             )
dataCol   <- c(
                     "GLASS" = "cyan3", 
                     "IVYGAP" = "chartreuse1", 
                     "HEILAND" = "magenta"
             )
corMethod <- "spearman"
```

## Calculate pairwise correlations

  - In this section, we construct a Spatial Omics Set (`soSet`) by
  integrating the cohort-specific cNMF objects and the corresponding
  rank specifications.
  - Internally, `SOSet()` concatenates the metagene loading matrices (the
  *W* matrices) across cohorts and selected ranks, thereby placing all
  metagenes into a common representation suitable for cross-dataset
  comparison.
  - It then computes an all-pairs correlation matrix across the
  concatenated metagenes using the specified correlation metric
  (`corMet`; here, Spearman).
  - The resulting `soSet` object stores the inputs, cohort metadata
  (including visualization colors), and the computed correlation
  structure, which serves as the basis for downstream network inference
  and community detection.

``` r
soSet <- SOSet(
        NMFobjL = dataL, 
        NMFrankL = rankL, 
        dataCol = dataCol, 
        corMet = corMethod
)
```

    ## 3 dataset(s) found in the list: GLASS, IVYGAP, HEILAND

    ## Loading: GLASS [Rank: 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
    ## Loading: IVYGAP [Rank: 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
    ## Loading: HEILAND [Rank: 5, 10, 15, 20, 25, 30]

    ## WARNING::Number of genes are different across datasets.

    ## Calculating: all pairwise correlation coefficients.

    ## Correlation computed with
    ## • Method: 'spearman'
    ## • Missing treated using: 'pairwise.complete.obs'
    ## 
    ## Assigned color(s):
    ##  > GLASS: cyan3
    ##  > IVYGAP: chartreuse1
    ##  > HEILAND: magenta

``` r
soSet
```

    ## Dataset(s): 
    ##   GLASS (cyan3)
    ##   IVYGAP (chartreuse1)
    ##   HEILAND (magenta)
    ## Select rank(s): 
    ##   GLASS: 3, 4, 5, 6, 7, 8, 9, 10, 15, 20
    ##   IVYGAP: 3, 4, 5, 6, 7, 8, 9, 10, 15, 20
    ##   HEILAND: 5, 10, 15, 20, 25, 30
    ## Basis (W) matrices: 
    ##   #Genes     : 26316
    ##   #Metagenes : 279
    ## Correlation method: 
    ##   spearman (pairwise.complete.obs)
    ## Correlation matrix: 
    ##   Symmetric matrix with 279 columns X 279 rows.

## Generate a correlation network

  - This section converts the correlation structure stored in `soSet` into
  a metagene similarity network and performs community detection under
  user-defined settings.
    1. First, a correlation threshold (`coefThre`) is applied to retain only
  sufficiently strong metagene-metagene associations (here,
  correlations ≥ 0.3), yielding a sparse graph representation.
    2. Community detection is then run with a fixed random seed (`seed`) to
  ensure reproducibility and a specified number of iterations (`niter`)
  to stabilize the optimization.
    3. The parameters `commWeight` and `cohortWeight` are _not_ used during 
    the community detection step. Instead, they are applied after communities 
    have been identified to rewire or reweight edges for visualization, 
    thereby influencing the spatial arrangement of nodes in network layouts.
    Specifically, `cohortWeight` increases the tendency for metagenes from 
    the same cohort (dataset) to cluster together, which is useful when 
    cohort-specific structure is a primary interpretive focus. In contrast, 
    `commWeight` increases the tendency for metagenes assigned to the same 
    community to cluster, which is preferable when emphasizing community-level 
    modular organization. By tuning these weights, users can generate 
    complementary visual representations of the same inferred communities 
    without altering the underlying community assignments.
    4. The resulting `soObj` stores the inferred network and the identified
  communities for downstream annotation and visualization.

``` r
corrCoefThre <- 0.3
seed         <- 1234
niter        <- 1000
commWeight   <- 100
cohortWeight <- 10

soObj <- SOTK(
        SOSet = soSet, 
        coefThre = corrCoefThre, 
        seed = seed, 
        niter = niter, 
        commWeight = commWeight, 
        cohortWeight = cohortWeight
)
```

    ## Seed: 1234

    ## Community search algorithm:
    ##  Fast Greedy

    ## Updating weights for
    ##  Community #1
    ##  Community #2
    ##  Community #3
    ##  Community #4
    ##  Community #5
    ##  Community #6
    ##  Community #7
    ##  Community #8

    ## Updating weights for
    ##  Data: GLASS
    ##  Data: IVYGAP
    ##  Data: HEILAND

    ## Calculating new layout based on new weights.

    ## Community-level network generated.

``` r
soObj
```

    ## Correlation network:
    ##   Nodes       : 279
    ##   Communities : 8 identified
    ## Parameters:
    ##   coefThre    : 0.3
    ##   seed        : 1234
    ##   niter       : 1000
    ##   drop        : FALSE
    ##   searchMet   : greedy
    ##   commWeight  : 100
    ##   cohortWeight: 10

## Save the object for reuse

  - To avoid repeating computationally intensive steps in later parts of
  the demo, we serialize the resulting `soObj` object to disk as an
  `.RDS` file.
  - This saved object can be reloaded in downstream sections to reproduce
  the same network and community assignments without rerunning the
  correlation-network construction and community detection.

``` r
saveRDS(soObj, file.path(download_dir, "soObj.RDS"))
message("soObj.RDS was created.")
```

    ## soObj.RDS was created.
