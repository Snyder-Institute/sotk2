#' Import consensus NMF (cNMF) results and convert to an NMF object
#'
#' Reads cNMF output files (gene spectra and consensus usages) across a range of factorization ranks
#' and converts them into an \code{"NMF.rank"} object compatible with downstream \pkg{NMF} workflows.
#'
#' The function expects the standard cNMF output file naming convention:
#' \itemize{
#'   \item \code{<prefix>.gene_spectra_tpm.k_<k>.dt_<denThre>.txt}
#'   \item \code{<prefix>.gene_spectra_score.k_<k>.dt_<denThre>.txt}
#'   \item \code{<prefix>.usages.k_<k>.dt_<denThre>.consensus.txt}
#' }
#'
#' @param prefix \code{character}. The \code{run_name} parameter used in the cNMF run.
#' @param metrics \code{character}. Either \code{"score"} or \code{"TPM"} to select which gene spectra file
#'   is used to construct the basis matrix \code{W}. Default is \code{"score"}.
#' @param minRank \code{numeric}. Minimum rank (\code{k}) to import.
#' @param maxRank \code{numeric}. Maximum rank (\code{k}) to import.
#' @param by \code{numeric}. Increment for the sequence of ranks. Default is 1.
#' @param denThre \code{character}. Density threshold string used by cNMF in file names (e.g., \code{"2_0"}).
#' @param cnmfDir \code{character}. Directory containing cNMF result files.
#' @param verbose \code{logical}. If \code{TRUE}, prints progress and matrix dimensions. Default is \code{TRUE}.
#'
#' @return An object of class \code{"NMF.rank"} with \code{$fit} containing imported fits indexed by \code{k}.
#'
#' @examples
#' nmfObj <- importCNMF(prefix = "PREFIX", metrics = "score", minRank = 2, maxRank = 10, by = 1, denThre = "2_0", cnmfDir = "./cNMF")
#'
#' @export
#' @importFrom NMF nmfModel
#' @importFrom methods new
importCNMF <- function(prefix = NULL, metrics = "score",
                       minRank = NULL, maxRank = NULL, by = 1, denThre = NULL,
                       cnmfDir = "NULL", verbose = TRUE) {
        if (is.null(prefix) || is.null(minRank) || is.null(maxRank) || is.null(denThre) || is.null(cnmfDir)) {
                stop("ERROR::Please provide all required arguments: prefix, minRank, maxRank, denThre, cnmfDir")
        }

        if (!any(metrics %in% c("score", "TPM"))) {
                stop("ERROR::Please specify either score or TPM")
        }

        if (!is.numeric(minRank) || !is.numeric(maxRank) || !is.numeric(by) || by < 1) {
                stop("ERROR::Please provide integer for both minRank/maxRank/by and/or by should be greater than 1")
        }

        resObj <- list()
        hitK <- c()
        for (k in seq(minRank, maxRank, by)) {
                tpm_path <- file.path(cnmfDir, paste0(prefix, ".gene_spectra_tpm.k_", k, ".dt_", denThre, ".txt"))
                score_path <- file.path(cnmfDir, paste0(prefix, ".gene_spectra_score.k_", k, ".dt_", denThre, ".txt"))
                usage_path <- file.path(cnmfDir, paste0(prefix, ".usages.k_", k, ".dt_", denThre, ".consensus.txt"))

                if (file.exists(tpm_path) &&
                        file.exists(score_path) &&
                        file.exists(usage_path)) {

                        if (verbose) message(paste0("\nLoading cNMF result files at k = ", k))

                        if (metrics == "TPM") {
                                wUsage <- t(read.table(tpm_path, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE))
                        } else if (metrics == "score") {
                                wUsage <- t(read.table(score_path, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE))
                        }

                        colnames(wUsage) <- NULL
                        if (verbose) message(paste("\tUsage (W): ", nrow(wUsage), " genes X ", ncol(wUsage), sep = ""))

                        hScore <- t(read.table(usage_path, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE))
                        rownames(hScore) <- NULL
                        if (verbose) message(paste("\tScore (H): ", nrow(hScore), " X ", ncol(hScore), " samples", sep = ""))

                        fitRes <- methods::new("NMFfitX1")
                        attr(fitRes, "fit") <- NMF::nmfModel(model = "NMFstd", W = wUsage, H = hScore)

                        if (length(resObj) < 1) {
                                resObj <- list(fitRes)
                        } else {
                                resObj <- append(resObj, fitRes)
                        }
                        hitK <- c(hitK, k)
                } else {
                        warning(paste0("WARNING::cNMF result files (i.e., *gene_spectra_*, *usage*) do not exist at k = ", k))
                }
        }

        if (length(hitK) > 0) {
                names(resObj) <- as.character(hitK)

                cnmf <- list(fit = resObj)
                class(cnmf) <- "NMF.rank"
        } else {
                stop("ERROR::No result found.")
        }

        return(cnmf)
}

#' Merge multiple NMF objects into one
#'
#' Combines multiple NMF fit objects into a single \code{NMF.rank}-like object by:
#' \enumerate{
#'   \item Computing quality measures across the input fits using \code{NMF::compare()}.
#'   \item Computing consensus connectivity matrices for each fit using \code{NMF::connectivity()}.
#'   \item Storing the original fits under \code{$fit}.
#' }
#'
#' @param nmfObjL \code{list}. A non-empty list of NMF fit objects compatible with \code{NMF::compare()} and
#'   \code{NMF::connectivity()}.
#'
#' @return An object of class \code{"NMF.rank"} with components:
#' \describe{
#'   \item{measures}{A matrix/data.frame of comparison measures (subset of columns as implemented).}
#'   \item{consensus}{A list of connectivity matrices (one per input fit).}
#'   \item{fit}{The original \code{nmfObjL} list.}
#' }
#'
#' @examples
#' nmfObj <- mergeNMFObjs(nmfObjL = objL)
#'
#' @export
#' @importFrom NMF compare connectivity
mergeNMFObjs <- function(nmfObjL) {
        if (!is.list(nmfObjL) || length(nmfObjL) < 1) {
                stop("ERROR::'nmfObjL' must be a non-empty list of NMF objects.")
        }

        measures <- compare(nmfObjL)
        measures <- measures[, c(5:ncol(measures))]
        rownames(measures) <- c(1:nrow(measures))

        consensus <- lapply(nmfObjL, function(res) {
                cons <- connectivity(res)
                attr(cons, "model") <- NULL
                return(cons)
        })

        fit <- nmfObjL

        obj <- list(
                measures = measures,
                consensus = consensus,
                fit = fit
        )
        class(obj) <- "NMF.rank"

        return(obj)
}

#' Extract metagene-associated genes (MAGs) from a correlation-network object
#'
#' Identifies metagene-associated genes (MAGs) for each metagene using a two-step procedure:
#' \enumerate{
#'   \item Select the top genes ranked by the basis (W) weights for each metagene.
#'   \item For those genes, compute the Spearman correlation between gene expression (per sample)
#'         and the corresponding coefficient (H) usage values (per sample), and keep genes with
#'         correlation greater than a threshold.
#' }
#' Source: Tsukahara, T. et al. A transcriptional rheostat couples past activity to future sensory responses. Cell 184, 6326–6343.e32 (2021).
#'
#' @param object A \code{SOSet} object.
#' @param profiles A named list of expression matrices. Each matrix should have genes as rows,
#'   samples as columns, and row names corresponding to gene identifiers.
#' @param noGenes \code{numeric} Maximum number of top-ranked genes to evaluate as MAG candidates.
#'   Default is 200.
#' @param coefficient \code{numeric} Spearman correlation threshold. Genes with correlation strictly greater than
#'   this value are retained. Default is 0.2.
#'
#' @return A nested list of MAGs. The top level is \code{dataName}, then \code{community}, then \code{metagene},
#'   each containing a character vector of retained genes.
#'
#' @examples
#' expr <- readRDS("GLASS_exprProfile.RDS")
#' mags <- getMAGs(soObj, list(GLASS = expr))
#'
#' @export
#' @importFrom stringr str_split
#' @importFrom igraph V
getMAGs <- function(object, profiles, noGenes = 200, coefficient = 0.2) {
        if (!is.list(profiles) || is.null(names(profiles))) {
                stop("ERROR::'profiles' must be a named list of expression matrices.")
        }
        if (!is.numeric(noGenes) || length(noGenes) != 1L || is.na(noGenes) || noGenes <= 0) {
                stop("ERROR::'noGenes' must be a single positive numeric value.")
        }
        if (!is.numeric(coefficient) || length(coefficient) != 1L || is.na(coefficient) || coefficient < 0 || coefficient > 1) {
                stop("ERROR::'coefficient' must be a single numeric value between 0 and 1 (inclusive).")
        }

        dataNames <- names(object@SOSet@NMFobjL)
        coefficientAcrossGEPs <- list()
        for (dataName in dataNames) {
                message("\nProcessing: ", dataName)
                if (dataName %in% names(profiles)) {
                        profile <- profiles[[dataName]]
                        rownames(profile) <- gsub("-", "\\.", rownames(profile))
                        
                        nmfRes <- object@SOSet@NMFobjL[[dataName]]
                        nodes <- igraph::V(object@corNetwork)
                        community <- nodes$community
                        names(community) <- nodes$name
                        subComm <- community[which(sapply(stringr::str_split(names(community), "\\$"), "[[", 1) == dataName)]

                        # Signature genes from each metagenes
                        for (eachCommunity in c(1:max(community))) {
                                coefficientAcrossGEPs[[dataName]][[eachCommunity]] <- list()
                                subsetCom <- names(subComm)[which(subComm == eachCommunity)]

                                for (comm in subsetCom) {
                                        message(paste0("Processing: Community #", eachCommunity, " from ", dataName, ", Metagenes: ", comm))

                                        buff <- unlist(stringr::str_split(comm, "\\$"))
                                        rank <- as.numeric(buff[2])
                                        k <- as.numeric(buff[3])

                                        H <- NMF::coef(get(as.character(k), nmfRes$fit))
                                        usageH <- as.matrix(H[rank,])

                                        W <- NMF::basis(get(as.character(k), nmfRes$fit))
                                        usageW <- W[,rank]
                                        usageW <- usageW[order(usageW, decreasing=T)]
                                        
                                        topGenes <- names(usageW)[1:noGenes]

                                        coefficients <- c()
                                        for (gene in topGenes) {
                                                df <- data.frame(
                                                        Sample = colnames(profile),
                                                        Expression = as.numeric(profile[which(rownames(profile) == gene),])
                                                )
                                                M <- merge(df, usageH, by.x = "Sample", by.y = "row.names", all=F)
                                                coefficients <- c(coefficients, cor(M[,2], M[,3], method = "spearman"))
                                        }
                                        names(coefficients) <- topGenes
                                        
                                        coefficientAcrossGEPs[[dataName]][[eachCommunity]][[comm]] <- names(coefficients)[coefficients > coefficient]
                                }
                        }
                } else {
                        message("WARNING::", dataName, " expression profile was not found in profiles.")
                }
        }
        return(coefficientAcrossGEPs)
}

#' Select contributing community genes from MAGs
#'
#' Aggregates metagene-associated genes (MAGs) at the community level and selects genes that
#' repeatedly appear across metagenes in the same community.
#'
#' A gene is selected as a contributing community gene if it appears as a MAG in more than
#' \code{proportion} of metagenes within that community.
#'
#' @param object A \code{SOSet} object
#' @param magL Nested list of MAGs returned by \code{getMAGs()}.
#' @param proportion Numeric. Threshold between 0 and 1 (inclusive). A gene must appear in
#'   more than this fraction of metagenes within a community to be retained. Default is 0.5.
#'
#' @return Named list where each element corresponds to a community (e.g., \code{"Community_1"}).
#'   Each community element is a named list keyed by dataset name (e.g., \code{dataName}) and
#'   contains a character vector of contributing genes.
#' @examples
#' mags <- getMAGs(soObj, list(GLASS = expr))
#' commGenes <- contributingCommunityGenes(soObj, mags)
#' 
#' @export
#' @importFrom igraph V
contributingCommunityGenes <- function(object, magL, proportion = 0.5) {

        if (!is.numeric(proportion) || length(proportion) != 1L || is.na(proportion) || proportion < 0 || proportion > 1) {
                stop("ERROR::'proportion' must be a single numeric value between 0 and 1 (inclusive).")
        }
        message("Threhold of the intersection was set to ", proportion)

        dataNames <- names(object@SOSet@NMFobjL)
        nodes <- igraph::V(object@corNetwork)
        community <- nodes$community

        if (is.null(community)) {
                stop("ERROR::Vertex attribute 'community' was not found in object@corNetwork.")
        }
        if (length(community) == 0) {
                stop("ERROR::Vertex attribute 'community' is empty in object@corNetwork.")
        }
        if (all(is.na(community))) {
                stop("ERROR::Vertex attribute 'community' contains only NA values in object@corNetwork.")
        }

        commSignatureGenes <- list()
        for (dataName in dataNames) {
                if (dataName %in% names(magL)) {
                        for (eachCommunity in c(1:max(community, na.rm = TRUE))) {
                                commSignatureGenes[[dataName]][[eachCommunity]] <- list()
                                if (length(magL[[dataName]][[eachCommunity]]) == 1) {
                                        genes <- magL[[dataName]][[eachCommunity]]
                                        outMat <- as.matrix(genes[[1]])
                                        outMat[,1] <- 1
                                } else {
                                        outMat <- as.matrix(table(unlist(magL[[dataName]][[eachCommunity]])))
                                        outMat <- as.matrix(outMat[order(outMat[,1], decreasing=T),])
                                }
                                if (length(which(outMat[,1]/length(names(magL[[dataName]][[eachCommunity]])) > proportion)) > 0) {
                                        commSignatureGenes[[dataName]][[eachCommunity]] <- rownames(outMat)[which(outMat[,1]/length(names(magL[[dataName]][[eachCommunity]])) > proportion)]
                                } else {
                                        message("\nWARNING::No intersected genes across metagenes greater than ", proportion, " in Community #", eachCommunity)
                                }
                        }
                } else {
                        message("\nWARNING::MAGs for ", dataName, " was not found in magL.")
                }
        }
        return(commSignatureGenes)
}

#' Get the most variable genes (MVGs)
#'
#' Selects the most variable genes from an expression matrix using a two-step filter:
#' \enumerate{
#'   \item Remove genes with coefficient of variation (CV) below \code{coefVar}.
#'   \item Rank remaining genes by standard deviation across samples and return the top \code{no} genes.
#' }
#'
#' @param profile A \code{matrix} of expression values where rows are genes and columns are samples.
#'   Row names must be gene identifiers.
#' @param coefVar \code{numeric} Coefficient of variation threshold. Genes with CV strictly lower than
#'   this value are removed. Default is 0.1.
#' @param no \code{numeric} Number of genes to return after ranking by standard deviation. Default is 2000.
#'
#' @return A character vector of gene identifiers (row names) corresponding to the selected MVGs.
#'
#' @examples
#' mvgs <- getMVGs(profile = as.matrix(expr), coefVar = 0.1, no = 2000)
#' # create a heatmap with MVGs
#'
#' @export
#' @importFrom EnvStats cv
getMVGs <- function(profile = NULL, coefVar = 0.1, no = 2000) {

        if (is.null(profile)) {
                stop("ERROR::Please provide a expression profile where rows are genes and columns are samples")
        }

        if (!is.matrix(profile)) {
                stop("ERROR::'profile' must be a matrix where rows are genes and columns are samples")
        }
        if (is.null(rownames(profile))) {
                stop("ERROR::'profile' must have row names corresponding to gene identifiers")
        }

        if (!is.numeric(coefVar) || length(coefVar) != 1L || is.na(coefVar) || coefVar < 0) {
                stop("ERROR::Please provide a non-negative numeric value for coefVar")
        }
        if (!is.numeric(no) || length(no) != 1L || is.na(no) || no <= 0) {
                stop("ERROR::Please provide a positive numeric value for no")
        }

        coefficientVariation <- apply(profile, 1, EnvStats::cv)
        if (length(which(coefficientVariation < coefVar)) > 0) {
                profile <- profile[-which(coefficientVariation < coefVar), ]
        }
        if (nrow(profile) < 1) {
                stop("ERROR::Please adjust coefVar")
        }

        standardDeviation <- apply(profile, 1, stats::sd)
        profile <- profile[order(standardDeviation, decreasing=T),]

        nOut <- min(as.integer(no), nrow(profile))
        return(rownames(profile)[seq_len(nOut)])
}
