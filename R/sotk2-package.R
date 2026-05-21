#' sotk2: cross-platform omics integration through deconvolution-derived modules
#'
#' @description
#' `sotk2` integrates NMF/cNMF outputs across datasets and modalities by
#' concatenating basis (W) matrices, computing metagene-metagene correlations,
#' and applying community detection on the resulting graph.
#'
#' @keywords internal
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics hist legend lines par plot.new title
#' @importFrom stats aggregate chisq.test cor density sd
"_PACKAGE"

# Declare ggplot aesthetic names as known globals.
utils::globalVariables(c("Community", "Proportion", "rank"))
