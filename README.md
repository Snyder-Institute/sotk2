![SOTK_logo](./.github/logo.png)

## About
### What is Spatial Omics Toolkit (SOTK)?
 * The Spatial Omics Toolkit (SOTK) is an R package that offers a comprehensive suite of functions for identifying biologically meaningful modules from spatial transcriptomics data. It enables integrative analysis across multiple cohorts or platforms by constructing a unified correlation network based on deconvoluted outputs, leveraging community detection algorithms to uncover overrepresented biological patterns.
 * SOTK also facilitates optimal rank selection from multiple non-negative matrix factorization (NMF) runs. Specifically, it identifies the minimal rank that captures the majority—if not all—of the communities identified via community search algorithms, providing a data-driven criterion for selecting latent factors.
 * Originally developed for Nanostring GeoMx Digital Spatial Profiler (DSP) data in 2021, SOTK is especially effective when multiple spatial segments (e.g., cell types) are independently profiled from the same tissue sample. This design allows users to incorporate and analyze numerous segment-level profiles per sample or patient.
 
### Features
 * Identifies biologically meaningful latent factors from deconvoluted spatial transcriptomics data
 * Data-driven rank selection across multiple NMF runs
 * Optimized for systematic exploration of specific cell types or tissue segments
 * Enables spatio-temporal analysis of gene sets
 * Supports integrative, cross-cohort spatial omics analysis

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
devtools::install_github("lootpiz/SOTK")
```

## How-to
### Workflow
 *	Unsupervised deconvolution
    -	Perform non-negative matrix factorization (NMF) across multiple ranks.  
 *	Identification of biological modules (i.e., sets of metagenes)
    1. Compute pairwise correlation coefficients of metagenes across all ranks.
    2. Construct a correlation network using only positively correlated metagenes.
    3. Apply the fast greedy community detection algorithm to identify modules within the network.
    4. Select the minimum rank that captures the largest number of communities (biological modules).
    5. Assign AOIs or samples to each metagene based on NMF coefficient values.

### Demonstration
  * SOTK is also available as an interactive ShinyApp for exploring GeoMx DSP data. You can try it here: [https://shinyapps.ucalgary.ca/SOTK](https://shinyapps.ucalgary.ca/SOTK)

