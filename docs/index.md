# Spatial Omics Toolkit 2 (sotk2)

- Last updated: December 27, 2025

---

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

### Key differences between SOTK and sotk2

**sotk2** extends the original SOTK workflow from “rank selection within a dataset” to **cross-dataset, cross-platform module integration**. The emphasis shifts from choosing an optimal *k* in a single analysis to building a comparable module space across datasets and extracting communities that persist across ranks and cohorts.

| Feature | sotk2 | SOTK | 
|---|---|---|
| Primary goal | Cross-dataset and cross-modality integration using deconvolution-derived modules | Data-driven selection of an optimal latent factor number within a dataset | 
| Input data types | Any modality with NMF/cNMF outputs (bulk, single-cell, spatial transcriptomics, protein abundance) | Spatial transcriptomics–focused | 
| Integration scope | Across datasets, cohorts, platforms, and modalities via correlation networks | Within dataset (single platform) | 
| Network representation | Metagene-level networks plus **community-level abstraction** to scale integration and interpretation | Metagene-level networks | 
| Composition assessment | **Residual-based overrepresentation** to summarize sample-type composition at the community level | Limited or dataset-specific | 
| Comparative visualization | Consistent-layout community networks to support direct cross-dataset comparison | Basic plotting | 
| Gene interpretation | Built-in extraction of **metagene-associated genes (MAGs)** and selection of **contributing community genes** for annotation | Metagene genes (limited) | 
| Intended outcome | Identify robust communities/modules and compare their presence and composition across datasets | Choose *k* and interpret metagenes | 

## Citation
If you use **sotk2** in your work, please cite the associated manuscript: _to be added_
