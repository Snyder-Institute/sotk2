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

### Key differences between SOTK and sotk2

**sotk2** extends the original SOTK workflow from “rank selection within a dataset” to **cross-dataset, cross-platform module integration**. The emphasis shifts from choosing an optimal *k* in a single analysis to building a comparable module space across datasets and extracting communities that persist across ranks and cohorts.

| Feature | SOTK | sotk2 |
|---|---|---|
| Primary goal | Data-driven selection of an optimal latent factor number within a dataset | Cross-dataset and cross-modality integration using deconvolution-derived modules |
| Input data types | Spatial transcriptomics–focused | Any modality with NMF/cNMF outputs (bulk, single-cell, spatial transcriptomics, protein abundance) |
| Integration scope | Within dataset (single platform) | Across datasets, cohorts, platforms, and modalities via correlation networks |
| Network representation | Metagene-level networks | Metagene-level networks plus **community-level abstraction** to scale integration and interpretation |
| Composition assessment | Limited or dataset-specific | **Residual-based overrepresentation** to summarize sample-type composition at the community level |
| Comparative visualization | Basic plotting | Consistent-layout community networks to support direct cross-dataset comparison |
| Gene interpretation | Metagene genes (limited) | Built-in extraction of **metagene-associated genes (MAGs)** and selection of **contributing community genes** for annotation |
| Intended outcome | Choose *k* and interpret metagenes | Identify robust communities/modules and compare their presence and composition across datasets |

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

Core imports typically include:
* igraph
* methods
* NMF
* RColorBrewer
* stringr

Additional packages may be used for plotting and summaries (for example, ggplot2, reshape2, scales) depending on which functions you call.

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

Below is a minimal example showing the object flow. For a full end-to-end script, see `inst/scripts/` in this repository and [Zenodo](https://doi.org/10.5281/zenodo.18063318).

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
The vignette is available online and provides step-by-step instructions for downloading the demo data, running sotk2, generating visualizations, and performing annotations to support interpretation. Please visit: <a href="https://Snyder-Institute.github.io/sotk2/" target="_blank">https://Snyder-Institute.github.io/sotk2/</a>


## Project provenance

My access to the SOTK repository was revoked. As a result, I did not fork the current SOTK GitHub repository; instead, I developed **sotk2** as an independent codebase that extends the original concepts with additional functionality intended to benefit the broader community.

## Citation
If you use **sotk2** in your work, please cite the associated manuscript:
