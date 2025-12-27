# Introduction

  - This section outlines a practical workflow for extracting
  metagene-associated genes (MAGs) and prioritizing highly contributing
  genes for each inferred community.
  - Communities are defined on the metagene correlation network, and
  community-level gene importance is derived by aggregating information
  across the metagenes assigned to a given community.
  - Identification of MAGs is essential for biological interpretation,
  because these genes provide a direct link between latent metagene
  structure and measurable transcriptional programs.
  - Community-level MAG summaries can be compared against established
  marker sets and published gene signatures to support functional
  annotation, assess concordance with known biological states, and
  reduce the risk of over-interpreting purely network-driven community
  structure.

## Directory settings

  - This block defines the working directory used throughout the demo to
  store downloaded inputs and generated outputs.

``` r
# If you want to use a user-defined output directory,
# uncomment and set the download_dir parameter.

# download_dir <- "/path/to/download"

if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && 
        nzchar(download_dir)) {
        download_dir <- download_dir
} else {
        download_dir <- tools::R_user_dir("sotk2", "data")
}
```

## Load the spatial omics object

  - We load the previously generated `soObj` object, which contains the
  correlation network and community detection results produced in the
  earlier steps of the workflow.
  - The script first attaches the `sotk2` package and then checks whether
  `soObj.RDS` is present in `download_dir`.
  - We then import the GLASS cohort expression profile (`expr_GLASS.RDS`) downloaded from Zenodo, which is required for downstream extraction of metagene-associated genes (MAGs) and community-level gene summaries.

``` r
library(sotk2)

if (file.exists(file.path(download_dir, "soObj.RDS"))) {
        soObj <- readRDS(file.path(download_dir, "soObj.RDS"))
} else {
        stop("ERROR: the soObj.RDS file not found.")
}

if (file.exists(file.path(download_dir, "expr_GLASS.RDS"))) {
        expr <- readRDS(file.path(download_dir, "expr_GLASS.RDS"))
} else {
        stop("ERROR: the soObj.RDS file not found.")
}
```

## Get metagene-associated genes

- This section extracts metagene-associated genes (MAGs) by integrating
  the GLASS expression matrix with the metagene and community structure
  stored in `soObj`.
  1.  The function `getMAGs()` is then applied to compute, for each
      metagene, the set of genes most strongly associated with that
      metagene within the specified cohort.
  2.  The resulting object (`mags`) is organized hierarchically by
      cohort, community, and metagene, enabling users to inspect MAGs at
      multiple resolutions (cohort-wide, community-level, and individual
      metagene-level).
  3.  The final lines provide an example of accessing MAGs for a
      specific metagene within a given community.

``` r
mags <- getMAGs(soObj, list(GLASS = expr))
message("MAGs for the metagenes in Community #1: ", names(mags[["GLASS"]][[1]])[1])
```

    ## MAGs for the metagenes in Community #1: GLASS$05$05

``` r
mags[["GLASS"]][[1]][[1]] # [[data/cohort]][[community #]][[metagenes]]
```

    ##   [1] "ADORA3"   "ADPRH"    "AGTRAP"   "AIF1"     "AKR1A1"   "ALDH3B1" 
    ##   [7] "ALOX5"    "ALOX5AP"  "ANXA2"    "APOBR"    "ARHGAP30" "ARHGDIB" 
    ##  [13] "ARPC1B"   "ARRB2"    "ASGR2"    "BATF"     "BRI3"     "C1QA"    
    ##  [19] "C1QB"     "C1QC"     "C1R"      "C1S"      "C2"       "C3AR1"   
    ##  [25] "CALHM2"   "CAP1"     "CAPG"     "CASP1"    "CASP4"    "CCR1"    
    ##  [31] "CCRL2"    "CD14"     "CD300A"   "CD300C"   "CD300LF"  "CD33"    
    ##  [37] "CD37"     "CD4"      "CD53"     "CD68"     "CD74"     "CEBPB"   
    ##  [43] "CEBPD"    "CERS2"    "CIB1"     "CLDN7"    "CLIC1"    "CNPY3"   
    ##  [49] "CSTB"     "CTSA"     "CTSB"     "CTSC"     "CTSD"     "CTSL"    
    ##  [55] "CTSZ"     "CXCL16"   "CXCR4"    "CYBA"     "DENND1C"  "DENND2D" 
    ##  [61] "DNASE2"   "DOK1"     "DOK2"     "DOK3"     "DUSP1"    "DUSP23"  
    ##  [67] "EDEM2"    "EFEMP2"   "EMP3"     "ENG"      "EVA1B"    "FAH"     
    ##  [73] "FBP1"     "FCER1G"   "FCGR2B"   "FCGRT"    "FERMT3"   "FLII"    
    ##  [79] "FTH1"     "FTL"      "FUCA2"    "GAA"      "GNA15"    "GPSM3"   
    ##  [85] "GPX1"     "GRN"      "GYPC"     "HCK"      "HCST"     "HEXA"    
    ##  [91] "HLA.B"    "HLA.DMA"  "HLA.DPA1" "HLA.DPB1" "HLA.DRA"  "HLA.DRB1"
    ##  [97] "HLA.E"    "HLX"      "HSD3B7"   "ICAM1"    "IER3"     "IFI30"   
    ## [103] "IFNGR2"   "IL10RB"   "IL15RA"   "ISG20"    "ITGB2"    "ITPRIP"  
    ## [109] "JUN"      "LAIR1"    "LAPTM5"   "LGALS1"   "LGALS3BP" "LILRB3"  
    ## [115] "LITAF"    "LRRC25"   "MAFB"     "MGAT1"    "MS4A4A"   "MS4A6A"  
    ## [121] "MSN"      "MVP"      "MXRA8"    "MYL9"     "MYO1G"    "NANS"    
    ## [127] "NCF2"     "NCF4"     "NFKBIA"   "NPC2"     "NUPR1"    "OSCAR"   
    ## [133] "PARVG"    "PFN1"     "PILRA"    "PLAUR"    "PLIN3"    "PLTP"    
    ## [139] "POLD4"    "PPM1M"    "PPP1R18"  "PTPN6"    "PYCARD"   "RAB20"   
    ## [145] "RAB42"    "RAC2"     "RASAL3"   "RASSF1"   "RCN3"     "RGS19"   
    ## [151] "RHOG"     "RILPL2"   "RIPK3"    "RNASE6"   "RNASET2"  "RNF130"  
    ## [157] "RRAS"     "S100A11"  "SASH3"    "SERF2"    "SERPINB1" "SIGLEC7" 
    ## [163] "SIGLEC9"  "SLAMF8"   "SLC10A3"  "SLC11A1"  "SLC15A3"  "SLC17A9" 
    ## [169] "SLC39A1"  "SPI1"     "SRGN"     "STAB1"    "STXBP2"   "SYNGR2"  
    ## [175] "TAGLN2"   "TAX1BP3"  "TCIRG1"   "TGFB1"    "THBD"     "THEMIS2" 
    ## [181] "TLN1"     "TMED9"    "TMEM150A" "TNFRSF14" "TNFRSF1B" "TNFSF13" 
    ## [187] "TNIP2"    "TRADD"    "TRIM21"   "TSPO"     "TWF2"     "TXNDC5"  
    ## [193] "TYMP"     "TYROBP"   "VAMP5"    "VAMP8"    "VAT1"     "VKORC1"  
    ## [199] "VSIG4"    "WAS"

## Get highly contributing genes

  - This section prioritizes community-level genes by aggregating
  metagene-associated gene information across all metagenes assigned to
  each community.
  - Using the MAGs extracted in the previous step,
  `contributingCommunityGenes()` identifies genes that contribute most
  strongly and consistently to a given community's metagene composition
  by quantifying how frequently each gene appears among the MAG sets of
  community member metagenes.
  - By default, the function applies a proportion threshold of `> 0.5`,
  meaning that a gene is retained as a highly contributing community
  gene only if it is identified as a MAG for more than half of the
  metagenes within that community.
  - For example, if a community contains 10 metagenes, a gene must be
  present among the MAGs of at least 6 metagenes to be selected under
  the default setting.
  - The resulting object is indexed by cohort and community, enabling
  direct inspection of top-ranked genes for each inferred module. The
  final lines print an example output showing the highly contributing
  genes for Community 1 in the GLASS cohort.

``` r
highlyContributingCommunityGenes <- contributingCommunityGenes(soObj, mags)
```

    ## Threhold of the intersection was set to 0.5

``` r
message("Highly contributing genes in Community #1 from GLASS:")
```

    ## Highly contributing genes in Community #1 from GLASS:

``` r
highlyContributingCommunityGenes[["GLASS"]][[1]] # Commuinty #1
```

    ##  [1] "ALOX5AP"  "ANXA2"    "C1R"      "C1S"      "CD14"     "CD300A"  
    ##  [7] "CD300C"   "CEBPB"    "CLIC1"    "CSTB"     "CTSB"     "CTSC"    
    ## [13] "CTSL"     "CTSZ"     "DOK2"     "FBP1"     "FCER1G"   "FCGRT"   
    ## [19] "GNA15"    "HCST"     "IER3"     "IFI30"    "LAPTM5"   "LILRB3"  
    ## [25] "MXRA8"    "MYO1G"    "NCF4"     "NFKBIA"   "OSCAR"    "PLAUR"   
    ## [31] "PLTP"     "RCN3"     "RRAS"     "S100A11"  "SLAMF8"   "STAB1"   
    ## [37] "TCIRG1"   "THBD"     "TNFRSF1B"


## Practical notes

### Methodological considerations: NMF, cNMF, and integration

End-to-end scripts—from NMF/cNMF execution to downstream `sotk2` analyses—are provided under `inst/scripts/`, enabling users to reproduce the full workflow and adapt it to their own datasets. Users may generate inputs using either standard NMF (https://github.com/renozao/NMF) or consensus NMF (cNMF; https://github.com/dylkot/cNMF). cNMF was originally developed for single-cell RNA-seq, where expression matrices are typically sparse. In this demonstration, we also apply cNMF to bulk RNA-seq because cNMF operates on the full expression matrix and returns factorization results over all genes. In contrast, standard NMF is often applied to a reduced feature space (for example, highly variable genes) to improve computational efficiency and to focus the factorization on dominant biological signals.

When a targeted gene panel is available (for example, molecular subtype signatures such as Verhaak and/or cell-state markers such as Neftel), bulk RNA-seq profiles may be subset accordingly prior to standard NMF, while cNMF can be reserved for sparse modalities such as single-cell or spatial transcriptomics. The resulting NMF and cNMF outputs can then be integrated in `sotk2` using correlation-based network construction and community detection.

### Analytical considerations: NMF/cNMF inputs and integration strategy

In this demonstration, ranks 3–10, 15, and 20 were evaluated for GLASS and IVYGAP (bulk RNA-seq), whereas spot-resolution Visium v1 data were evaluated at ranks 5–30 in increments of 5 (that is, 5, 10, 15, 20, 25, 30). In our experience, inclusion of very high ranks often has limited impact on the dominant community structure identified by `sotk2`, because higher-rank components can capture increasingly sample-specific variation, minor biological features, or outlier-driven structure. The optimal rank range is dataset-dependent and is influenced by sample size and heterogeneity.

As a practical starting point, we recommend prioritizing lower ranks while avoiding rank 2, which can be disproportionately sensitive to single outliers. Users should assess potential outliers by inspecting metagene usage profiles (and related diagnostic plots) before finalizing the rank grid. Notably, including higher ranks increases the number of metagenes available for network construction; this can be useful when the goal is to interrogate rare programs or outlier-associated structure, but it may also introduce repeated small communities that complicate interpretation.

### Network construction and integration parameters
The correlation threshold used to define edges in the integrated network is a key analytic parameter and should be selected empirically through exploratory evaluation (for example, by examining coefficient distributions and the resulting network sparsity and connectivity). Joint integration of NMF and cNMF outputs is supported; however, the effective shared gene set may be reduced when correlations are computed only over genes present in all inputs. Running NMF on the full expression matrix can lessen this constraint, although the associated computational burden may increase substantially with matrix size and the rank grid under consideration.

By default, network construction retains positive correlations exceeding a user-specified cutoff. If anti-correlated structure is scientifically relevant, an alternative thresholding scheme should be defined a priori and explicitly documented (for example, thresholding on absolute correlation while preserving edge sign). The “weighted layout” option is a post hoc visualization procedure that re-weights or rewires edges to encourage co-localization of nodes from the same community and/or cohort, improving graphical interpretability. Importantly, this layout strategy does not alter the underlying correlation estimates or the community assignments.