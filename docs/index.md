# Spatial Omics Toolkit 2 (sotk2)

## What is sotk2
**sotk2** is an R package for integrating omics datasets using modules derived from non-negative matrix factorization (NMF) or consensus NMF (cNMF). The core idea is to treat each gene expression program (metagene) as a comparable unit across datasets, then integrate programs through a correlation-based network followed by community detection.

**sotk2** is designed to be self-contained and independent: it does not require any prior packages or objects outside this repository. Inputs may come from any platform or modality (for example, bulk RNA-seq, single-cell RNA-seq, spatial transcriptomics, or protein abundance), as long as NMF/cNMF outputs are available (or can be imported).

### Features
 * Identification of biologically meaningful latent factors from deconvoluted omics data
 * Data-driven rank selection across multiple NMF or cNMF runs
 * Correlation-based integration of biological modules across datasets and platforms
 * Support for spatial, bulk, single-cell, and protein-level omics data
 * Community abstraction and network-level visualization for large integrative analyses
 * Assessing sample-type enrichment using Pearson residuals (observed vs expected counts; Chi-squared framework)
 * Consistent network layouts for comparative, cross-dataset interpretation

## Concepts
**sotk2** organizes analysis into two primary objects:

  * **SpatialOmicsSet**
    * Stores per-dataset NMF results (as `NMF.rank` objects)
    * Concatenates basis (W) matrices across ranks/datasets
    * Computes the metagene–metagene correlation matrix

  * **MetageneCorrelationNetwork**
    * Thresholds correlations and builds a metagene graph
    * Detects communities and stores community membership
    * Computes layouts for plotting and creates an optional community-level aggregated network

## Try it locally
The fastest way to see sotk2 in action is the bundled quickstart vignette, which runs on the demo data shipped with the package and requires no external downloads:

```r
install.packages("devtools")
devtools::install_github("Snyder-Institute/sotk2", build_vignettes = TRUE)
vignette("sotk2 quickstart", package = "sotk2")
```

A self-hosted ShinyApp and a portable Docker image are in preparation; this page will be updated when they are live.

## Citation
If you use **sotk2** in your work, please cite the associated manuscript: _to be added_
