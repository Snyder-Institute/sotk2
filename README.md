# Spatial Omics Toolkit 2 (SOTK2)
> Cross-Platform Omics Integration Through Deconvolution-Derived Modules

![Logo](.github/imgs/sotk2_logo_black.png#gh-dark-mode-only)
![Logo](.github/imgs/sotk2_logo_white.png#gh-light-mode-only)

## About
### What Is Spatial Omics Toolkit 2 (SOTK2)?
  * Spatial Omics Toolkit 2 (SOTK2) is a major upgrade of the original Spatial Omics Toolkit (SOTK). While SOTK focused on selecting the optimal number of latent factors from spatial transcriptomics data in a data-driven manner, **SOTK2 extends this framework to enable integrative analysis across heterogeneous omics datasets** through biologically informed module integration.
  * SOTK2 provides a correlation-based strategy to integrate biological modules derived from non-negative matrix factorization (NMF) or consensus NMF (cNMF) across datasets, cohorts, platforms, and modalities. Inputs are no longer limited to spatial transcriptomics data and can include bulk RNA sequencing, single-cell RNA sequencing, spatial transcriptomics, and protein expression profiles.
  * Originally developed for Nanostring GeoMx Digital Spatial Profiler (DSP) data in 2021, SOTK2 retains strong support for spatially resolved data while enabling cross-platform and cross-cohort integration through unified module-level representations.

### Features
 * Identification of biologically meaningful latent factors from deconvoluted omics data
 * Data-driven rank selection across multiple NMF or cNMF runs
 * Correlation-based integration of biological modules across datasets and platforms
 * Support for spatial, bulk, single-cell, and protein-level omics data
 * Community abstraction and network-level visualization for large integrative analyses
 * Quantification of sample-type overrepresentation using residual-based statistics
 * Consistent network layouts for comparative, cross-dataset interpretation

### Key Differences Between SOTK and SOTK2

| Feature | SOTK | SOTK2 |
|-------|------|-------|
| Primary Goal | Data-driven selection of optimal latent factor number | Cross-dataset and cross-modality data integration |
| Input Data Types | Spatial transcriptomics | Spatial transcriptomics, bulk/single-cell RNA-seq, protein profiles |
| Biological Module Integration | Within dataset | Across datasets, platforms, and modalities |
| Network Abstraction | Individual metagenes | Community-level nodes with residual-based composition |
| Visualization Support | Basic | Community abstraction, edge aggregation, and comparative layouts |

## Installation
### Requirement
 * R > 4.3.0
 * Dependencies (alphabetical order):
   * corrr
   * grid
   * igraph
   * methods
   * NMF
   * RColorBrewer
   * stringr

### Install via GitHub
```
install.packages("devtools")
devtools::install_github("Snyder-Institute/sotk2")
```

## How-to
### Workflow
1. **Unsupervised Deconvolution**
   - Perform NMF or cNMF (consensus NMF) independently for each dataset/cohort across multiple ranks.

2. **Concatenation of Gene Expression Programs**
   - Concatenate _W_ matrices across ranks and datasets to form a unified feature space.

3. **Correlation-Based Network Construction**
   - Compute pairwise correlations among metagenes.
   - Retain positively correlated edges using dataset-specific thresholds to mitigate batch effects.

4. **Community Detection**
   - Apply community detection algorithms, such as fast greedy clustering, to identify biological modules.

5. **Optimal Rank Selection**
   - Select the minimum rank that captures the majority of detected biological communities.

6. **Community Abstraction and Visualization**
   - Abstract metagenes into community-level nodes.
   - Scale node size by the number of assigned gene expression programs.
   - Scale edge thickness by the number of inter-community connections.

7. **Sample Composition and Residual Analysis**
   - Assign samples to gene expression programs based on maximal usage.
   - Quantify overrepresentation using residuals computed as:
     
     ``` residual = (observed - expected) / sqrt(expected) ```
     
   - Replace negative residuals with zero.
   - Visualize residual-based sample composition within community nodes.

8. **Cross-Dataset Comparison**
   - Generate community networks using consistent layouts to enable direct comparison of biological patterns across datasets.

### Demonstration
  * SOTK is also available as an interactive ShinyApp for exploring GeoMx DSP data. You can try it here: [https://shinyapps.ucalgary.ca/SOTK](https://shinyapps.ucalgary.ca/SOTK)

