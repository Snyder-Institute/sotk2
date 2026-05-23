# Spatial Omics Toolkit 2 (sotk2)
> Cross-platform omics integration through deconvolution-derived modules

![Logo](.github/imgs/sotk2_logo_black.png#gh-dark-mode-only)
![Logo](.github/imgs/sotk2_logo_white.png#gh-light-mode-only)

## About
### What is sotk2
**sotk2** is an R package for integrating omics datasets using **modules derived from non-negative matrix factorization (NMF) or consensus NMF (cNMF)**. The core idea is to treat each gene expression program (metagene) as a comparable unit across datasets, then integrate programs through a **correlation-based network** followed by **community detection**.

**sotk2** is designed to be **self-contained and independent**: it does not require any prior packages or objects outside this repository. Inputs may come from any platform or modality (for example, bulk RNA-seq, single-cell RNA-seq, spatial transcriptomics, or protein abundance), as long as NMF/cNMF outputs are available (or can be imported).

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

## Interactive demo

The Shiny app is distributed as a portable Docker image. The image bundles the demo data, so a single `docker run` is sufficient to explore the package end-to-end without local R setup. It is multi-arch (`linux/amd64` + `linux/arm64`), so the right binary is selected automatically on Intel/AMD Linux, Apple Silicon, and AWS Graviton hosts.

```bash
docker pull thebiohub/sotk2:1.0.0
docker run --rm -p 11630:11630 thebiohub/sotk2:1.0.0
open http://localhost:11630
```

Image registry: <a href="https://hub.docker.com/r/thebiohub/sotk2" target="_blank">hub.docker.com/r/thebiohub/sotk2</a>. The companion Dockerfile and launcher script live at <a href="https://github.com/Snyder-Institute/sotk2-docker" target="_blank">Snyder-Institute/sotk2-docker</a>.

If you prefer not to use Docker, two alternatives:

1. **Bundled quickstart vignette.** Single cohort, no external downloads, runs on the demo data shipped with the R package:

   ```r
   install.packages("devtools")
   devtools::install_github("Snyder-Institute/sotk2", build_vignettes = TRUE)
   vignette("sotk2 quickstart", package = "sotk2")
   ```

2. **Online walkthrough.** The documentation at <a href="https://Snyder-Institute.github.io/sotk2/" target="_blank">https://Snyder-Institute.github.io/sotk2/</a> reproduces the GLASS / IVYGAP / HEILAND multi-cohort integration step by step, including all figures.

A hosted, self-managed ShinyApp deployment is in preparation; this section will be updated when it goes live.

---

## Installation

### Requirements
* R >= 4.3.0

### Install from GitHub
```r
install.packages("devtools")
devtools::install_github("Snyder-Institute/sotk2")
```

### Dependencies
The package uses standard R infrastructure plus several common analysis/visualization packages. Exact versions are tracked in `DESCRIPTION`.

Imports:
* corrr
* EnvStats
* ggplot2
* igraph
* methods
* NMF
* RColorBrewer
* reshape2
* scales
* stats
* stringr

## Workflow overview

1. **Run NMF or cNMF per dataset**
   * Run NMF/cNMF independently for each dataset (optionally across a range of ranks).

2. **Concatenate basis matrices across ranks and datasets**
   * Combine W matrices into a unified basis matrix where each column corresponds to a metagene.

3. **Compute metagene correlations**
   * Compute pairwise correlations among metagenes using the concatenated basis matrix.

4. **Build a metagene correlation network**
   * Threshold correlations to keep robust positive edges (thresholds can be tuned per use case).

5. **Detect communities**
   * Apply community detection to identify groups of related metagenes (modules).

6. **Summarize and visualize**
   * Visualize metagene networks and community-level networks
   * Generate community summaries (counts, rank composition, and community-level connectivity)

7. **Gene-level interpretation (optional)**
   * Extract metagene-associated genes (MAGs)
   * Select contributing community genes for functional annotation

### Quick start

Below is a minimal example showing the object flow. For a full end-to-end script, see `inst/scripts/` in this repository and <a href="https://doi.org/10.5281/zenodo.18063318" target="_blank">Zenodo</a>.

```r
library(sotk2)

# Example inputs: per-dataset NMF.rank objects (one per dataset)
nmf_rank_object_list <- list(
  dataset_a = nmf_rank_dataset_a,
  dataset_b = nmf_rank_dataset_b
)

included_rank_list <- list(
  dataset_a = 3:10,
  dataset_b = 3:10
)

dataset_color_map <- c(
  dataset_a = "cyan3",
  dataset_b = "chartreuse1"
)

soSet <- SOSet(
        NMFobjL = nmf_rank_object_list, 
        NMFrankL = included_rank_list, 
        dataCol = dataset_color_map
)

soObj <- SOTK(SOSet = soSet)
```

---

## Vignette
The vignette is available online and provides step-by-step instructions for downloading the demo data, running sotk2, generating visualizations, and performing annotations to support interpretation. 

Access the **sotk2** vignette here: <a href="https://Snyder-Institute.github.io/sotk2/" target="_blank">https://Snyder-Institute.github.io/sotk2/</a>

## Background

**sotk2** is an independent implementation that extends concepts from earlier deconvolution-based integration toolkits into a cross-dataset, cross-modality framework with community-level abstraction, residual-based composition assessment, and built-in gene-level interpretation.

## Citation
If you use **sotk2** in your work, please cite:
> Reproducible transcriptional modules define glioblastoma ecosystems across independent cohorts. *bioRxiv* (2026). 
> https://doi.org/10.64898/2026.05.20.726700
 
