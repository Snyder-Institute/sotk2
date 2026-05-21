# sotk2 1.0.0

First public release of `sotk2`.

## Highlights

* `SOSet()` — S4 class for cross-dataset NMF/cNMF inputs and the concatenated
  metagene–metagene correlation matrix.
* `SOTK()` — S4 class for the thresholded correlation network, community
  assignments (greedy / leiden / betweenness / randomwalk / eigen / louvain),
  and a community-level aggregated network with consistent layouts for
  cross-cohort comparison.
* `importCNMF()` — load cNMF gene-spectra and consensus-usage outputs into an
  `NMF.rank`-compatible object so cNMF results can flow directly into the
  same `SOSet()` constructor used for `NMF` outputs.
* `mergeNMFObjs()` — combine independently-fit NMF runs (e.g., one per rank
  on different machines) into a single `NMF.rank` object suitable for
  `SOSet()`. Column selection from `NMF::compare()` is name-based and survives
  upstream changes in the `NMF` package.
* `plotCorrDensity()`, `plotNetwork()`, `plotCommNetwork()`, `selectRank()`,
  `statComm()` — visualization methods for correlation distributions,
  metagene-level networks, community-level networks, rank coverage across the
  network, and per-community rank composition. All accept `filename = NULL`
  to draw on the active device.
* `plotCommNetwork()` with `vertexInfo` renders pie annotations encoding
  **positive Pearson residuals** from a chi-squared test, summarizing
  sample-type or annotation over-representation per community.
* `getMAGs()` and `contributingCommunityGenes()` — gene-level interpretation:
  extract metagene-associated genes per metagene, then aggregate to select
  genes that contribute consistently across a community's metagenes.
* `getMVGs()` — utility for selecting the most variable genes from an
  expression matrix.
* `geoMean()` — geometric-mean helper used to summarize per-sample community
  activity from normalized metagene usages.
* `download_demo_data()` — fetch the Zenodo demo bundle
  ([10.5281/zenodo.18063318](https://doi.org/10.5281/zenodo.18063318))
  from R with MD5 verification (`set = "core"` for ~120 MB or
  `set = "full"` for ~660 MB, adding the Visium files).

## Documentation

* Online documentation aligned with the installed API
  (<https://Snyder-Institute.github.io/sotk2/>).

## Infrastructure

* `testthat` 3rd-edition test suite under `tests/testthat/`.
* GitHub Actions `R-CMD-check` workflow across Ubuntu, macOS, and Windows
  on R-release.
* `R CMD check` passes with no errors, warnings, or notes.
