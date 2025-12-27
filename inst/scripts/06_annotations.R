#!/usr/bin/env Rscript
# Once communities have been identified, users can annotate each community using sample-level metadata. 
# The script below generates annotation plots (e.g., primary vs. recurrence) and residual-based summaries. 
# Users can also compute the geometric mean of metagene usage (from the NMF results) to identify communities with stronger or weaker activity and to reduce the risk of misinterpretation.

# [WARNING] 
# The code below depends on the specific data structure and should be adapted for each dataset. 
# For example, in IVYGAP, sample annotations can be inferred directly from the sample names embedded in the NMF inputs, 
# whereas GLASS sample identifiers are barcodes that do not encode metadata. 
# Therefore, the annot_GLASS.RDS file is required to visualize sample-level annotations for the GLASS dataset.

# If you want to use a user-defined output directory, uncomment and set the download_dir parameter.
# download_dir <- "/path/to/download" # where soObj.RDS is located

if (exists("download_dir") && is.character(download_dir) && length(download_dir) == 1 && nzchar(download_dir)) {
        download_dir <- download_dir
} else {
        download_dir <- tools::R_user_dir("sotk2", "data")
}

if (!dir.exists(download_dir)) {
        dir.create(download_dir, recursive = TRUE)
        message(download_dir, " created.")
}

library(sotk2)

if (file.exists(file.path(download_dir, "soObj.RDS"))) {
        soObj <- readRDS(file.path(download_dir, "soObj.RDS"))
} else {
        stop("Run 04_runSotk2.R or download soObj.RDS from Zenodo.")
}

# Pie colors
pieColors <- list(
        IVYGAP = c("CT" = "red", "IT" = "orange", "LE" = "gold", "MVP" = "darkslategray4", "PAN" = "burlywood4"),
        GLASS = list(
                "sampleType" = c("Primary" = "chartreuse1", "Recurrence" = "darkgreen"),
                "molecularSubtype" = c("Classical" = "deepskyblue2", "Mesenchymal" = "deeppink3", "Proneural" = "coral3")
        )
)

# User defined function for annotation
.getIVYGAPidx <- function(x) {
        ctIdx <- c(); itIdx <- c(); leIdx <- c(); mvpIdx <- c(); panIdx <- c()

        if (!is.null(x)) {
                region <- stringr::str_sub(sapply(stringr::str_split(x, "__"), "[[", 2), 0, 3)
                ctIdx <- which(region == "Cel")
                itIdx <- which(region == "Inf")
                leIdx <- which(region == "Lea")
                mvpIdx <- which(region == "Mic")
                panIdx <- which(region == "Pse")
        }

        return(list(CT = ctIdx, IT = itIdx, LE = leIdx, MVP = mvpIdx, PAN = panIdx)) 
}

.getGLASSidx <- function(x, db) {
        priIdx <- c(); recIdx <- c()
        claIdx <- c(); mesIdx <- c(); proIdx <- c()

        if (!is.null(x)) {
                sub <- db[which(rownames(db) %in% x), c("sample_type", "Subtype")]
                priIdx <- which(sub$sample_type == "Primary")
                recIdx <- which(sub$sample_type == "Recurrence")
                claIdx <- which(sub$Subtype == "Classical")
                mesIdx <- which(sub$Subtype == "Mesenchymal")
                proIdx <- which(sub$Subtype == "Proneural")
        }

        return(list(
                "sampleType" = list(Pri = priIdx, Rec = recIdx), 
                "molecularSubtype" = list(Cla = claIdx, Mes = mesIdx, Pro = proIdx)
        ))
}

geoMean <- function(x) exp(mean(log(x)))

corNetwork <- soObj@corNetwork
clusterMembership <- soObj@sample2metagene
community <- c(1:length(soObj@commCols))
commNetwork <- soObj@commNetwork

allGEPs <- data.frame(
        Data = sapply(stringr::str_split(igraph::V(corNetwork)$name, "\\$"), "[[", 1), 
        GEP = igraph::V(corNetwork)$name,
        Community = igraph::V(corNetwork)$community
)

# Community annotation - GLASS
if (file.exists(file.path(download_dir, "annot_GLASS.RDS"))) {
        glassAnnot <- readRDS(file.path(download_dir, "annot_GLASS.RDS"))
        message("annot_GLASS.RDS file imported.")
} else {
        stop("Please download annot_GLASS.RDS by running 03_download.R.")
}

dName <- "GLASS"
legend <- names(pieColors[[dName]][["sampleType"]])
legendCol <- pieColors[[dName]][["sampleType"]]

cl <- clusterMembership[[dName]]
subGEPs <- allGEPs[which(allGEPs$Data == dName),]

# init
vertexLabel <- rep("", length(community)); names(vertexLabel) <- community
vertexSize <- rep(0.1, length(community)); names(vertexSize) <- community
vertexPie <- rep_len(list(numeric(length(legend))), length(community)); names(vertexPie) <- paste0("Community_", community)

for (whichComm in community) {
message(whichComm)
        commName <- paste0("Community_", whichComm)
        commSpecificGEPs <- subGEPs[which(subGEPs$Community == whichComm),]
        if (nrow(commSpecificGEPs) > 0) {
                allSamples <- c()
                for (gep in commSpecificGEPs$GEP) {
                        allSamples <- c(allSamples, cl[[gep]])
                }
                noSamples <- .getGLASSidx(unique(allSamples), glassAnnot)
                noSamples <- sapply(noSamples[["sampleType"]], length)
                
                if (sum(noSamples) != 0) {
                        vertexPie[[commName]] <- noSamples
                        vertexSize[whichComm] <- igraph::V(commNetwork)$size[whichComm]
                        vertexLabel[whichComm] <- igraph::V(commNetwork)$name[whichComm]
                } 
        }
}

vertexInfo <- list(vertexLabel = vertexLabel, vertexSize = vertexSize, vertexPie = vertexPie, legend = legend, legendCol = legendCol)
plotCommNetwork(soObj, vertexInfo = vertexInfo, filename = file.path(download_dir, "12_GLASS_sampleType_community_network.pdf"))

legend <- names(pieColors[[dName]][["molecularSubtype"]]) # Verhaak
legendCol <- pieColors[[dName]][["molecularSubtype"]]

# init
vertexLabel <- rep("", length(community)); names(vertexLabel) <- community
vertexSize <- rep(0.1, length(community)); names(vertexSize) <- community
vertexPie <- rep_len(list(numeric(length(legend))), length(community)); names(vertexPie) <- paste0("Community_", community)

for (whichComm in community) {
        commName <- paste0("Community_", whichComm)
        commSpecificGEPs <- subGEPs[which(subGEPs$Community == whichComm),]
        if (nrow(commSpecificGEPs) > 0) {
                allSamples <- c()
                for (gep in commSpecificGEPs$GEP) {
                        allSamples <- c(allSamples, cl[[gep]])
                }
                noSamples <- .getGLASSidx(unique(allSamples), glassAnnot)
                noSamples <- sapply(noSamples[["molecularSubtype"]], length)
                
                if (sum(noSamples) != 0) {
                        vertexPie[[commName]] <- noSamples
                        vertexSize[whichComm] <- igraph::V(commNetwork)$size[whichComm]
                        vertexLabel[whichComm] <- igraph::V(commNetwork)$name[whichComm]
                } 
        }
}

vertexInfo <- list(vertexLabel = vertexLabel, vertexSize = vertexSize, vertexPie = vertexPie, legend = legend, legendCol = legendCol)
plotCommNetwork(soObj, vertexInfo = vertexInfo, filename = file.path(download_dir, "13_GLASS_molecularType_community_network.pdf"))

## Community annotation - IVYGAP
dName <- "IVYGAP"
legend <- c("Cellular_Tumor", "Infiltrating_Tumor", "Leading_Edge", "Microvascular_proliferation", "Pseudopalisading_cells_around_necrosis")
legendLbl <- names(pieColors[[dName]])
legendCol <- pieColors[[dName]]

cl <- clusterMembership[[dName]]
subGEPs <- allGEPs[which(allGEPs$Data == dName),]

# init
vertexLabel <- rep("", length(community)); names(vertexLabel) <- community
vertexSize <- rep(0.1, length(community)); names(vertexSize) <- community
vertexPie <- rep_len(list(numeric(length(legend))), length(community)); names(vertexPie) <- paste0("Community_", community)

for (whichComm in community) {
        commName <- paste0("Community_", whichComm)
        commSpecificGEPs <- subGEPs[which(subGEPs$Community == whichComm),]
        if (nrow(commSpecificGEPs) > 0) {
                allSamples <- c()
                for (gep in commSpecificGEPs$GEP) {
                        allSamples <- c(allSamples, cl[[gep]])
                }
                noSamples <- .getIVYGAPidx(unique(allSamples))
                noSamples <- sapply(noSamples, length)
                
                if (sum(noSamples) != 0) {
                        vertexPie[[commName]] <- noSamples
                        vertexSize[whichComm] <- igraph::V(commNetwork)$size[whichComm]
                        vertexLabel[whichComm] <- igraph::V(commNetwork)$name[whichComm]
                } 
        }
}
vertexInfo <- list(vertexLabel = vertexLabel, vertexSize = vertexSize, vertexPie = vertexPie, legend = legend, legendLbl = legendLbl, legendCol = legendCol)
plotCommNetwork(soObj, vertexInfo = vertexInfo, filename = file.path(download_dir, "14_IVYGAP_community_network.pdf"))

# Community-level NMF usages
glass <- readRDS(file.path(download_dir, "nmfRes_GLASS.RDS"))
ivygap <- readRDS(file.path(download_dir, "nmfRes_IVYGAP.RDS"))
dataL <- list(GLASS = glass, IVYGAP = ivygap)

dName <- "GLASS"
nmfObj <- dataL[[dName]]
nodes <- igraph::V(soObj@corNetwork)
community <- nodes$community
names(community) <- nodes$name

pdf(file.path(download_dir, paste0("16_", dName, "_GEP_usage_geomean.pdf")), width=6, height=6)
for (comm in c(1:length(soObj@commCols))) {
        metagenes <- names(community)[which(community == comm)] # get metagenes in a given community
        metagenes <- metagenes[which(sapply(stringr::str_split(metagenes, "\\$"), "[[", 1) == dName)] # get metagenes in a given cohort

        if (length(metagenes) > 0) {
                usageMat <- c()
                for (metagene in metagenes) {
                        buff <- unlist(stringr::str_split(metagene, "\\$"))
                        k <- as.numeric(buff[3])
                        rank <- as.numeric(buff[2])
                        
                        usage <- NMF::coef(get(as.character(k), nmfObj$fit))
                        
                        rowSum <- apply(usage, 1, sum)
                        nUsage <- usage/rowSum
                        usageMat <- rbind(usageMat, nUsage[rank,])
                }

                usageMat[usageMat == 0] <- 1E-08
                rep <- apply(usageMat, 2, geoMean)
                names(rep) <- colnames(usageMat)

                idx <- .getGLASSidx(names(rep), glassAnnot)

                priIdx <- idx[["sampleType"]][["Pri"]]
                recIdx <- idx[["sampleType"]][["Rec"]]
                claIdx <- idx[["molecularSubtype"]][["Cla"]]
                mesIdx <- idx[["molecularSubtype"]][["Mes"]]
                proIdx <- idx[["molecularSubtype"]][["Pro"]]

                pri <- rep[priIdx]; pri <- pri[order(pri)]; priF <- ecdf(rep[priIdx])
                rec <- rep[recIdx]; rec <- rec[order(rec)]; recF <- ecdf(rep[recIdx])
                cla <- rep[claIdx]; cla <- cla[order(cla)]; claF <- ecdf(rep[claIdx])
                mes <- rep[mesIdx]; mes <- mes[order(mes)]; mesF <- ecdf(rep[mesIdx])
                pro <- rep[proIdx]; pro <- pro[order(pro)]; proF <- ecdf(rep[proIdx])
                
                plot(pri, priF(pri), col="chartreuse1", xlab="GeoMean(GEPs)", ylab="ECDF", main=paste0("Community #", comm), type="l", lwd=3, xlim=c(0, max(rep)), ylim=c(0, 1))
                lines(rec, recF(rec), col="darkgreen", lwd=3)
                lines(cla, claF(cla), col="deepskyblue2", lwd=3, lty=2)
                lines(mes, mesF(mes), col="deeppink3", lwd=3, lty=2)
                lines(pro, proF(pro), col="coral3", lwd=3, lty=2)
                
                legend("bottomright", 
                        legend=c(
                                paste0("Pri (", length(priIdx), ")"), 
                                paste0("Rec (", length(recIdx), ")"),
                                NA, 
                                paste0("Classical (", length(claIdx), ")"),
                                paste0("Mesenchymal (", length(mesIdx), ")"),
                                paste0("Proneural (", length(proIdx), ")")
                        ),
                        col = c("chartreuse1", "darkgreen", NA, "deepskyblue2", "deeppink3", "coral3"), 
                        pch = 15, pt.cex = 2.8,
                        bty="n", ncol=2
                )
        } else {
                message(paste0("No metagenes in Community #", comm))
        }
}
dev.off()

dName <- "IVYGAP"
nmfObj <- dataL[[dName]]
nodes <- igraph::V(soObj@corNetwork)
community <- nodes$community
names(community) <- nodes$name

pdf(file.path(download_dir, paste0("17_", dName, "_GEP_usage_geomean.pdf")), width=6, height=6)
for (comm in c(1:length(soObj@commCols))) {
        metagenes <- names(community)[which(community == comm)] # get metagenes in a given community
        metagenes <- metagenes[which(sapply(stringr::str_split(metagenes, "\\$"), "[[", 1) == dName)] # get metagenes in a given cohort

        if (length(metagenes) > 0) {
                usageMat <- c()
                for (metagene in metagenes) {
                        buff <- unlist(stringr::str_split(metagene, "\\$"))
                        k <- as.numeric(buff[3])
                        rank <- as.numeric(buff[2])
                        
                        usage <- NMF::coef(get(as.character(k), nmfObj$fit))
                        
                        rowSum <- apply(usage, 1, sum)
                        nUsage <- usage/rowSum
                        usageMat <- rbind(usageMat, nUsage[rank,])
                }
        
                usageMat[usageMat == 0] <- 1E-08
                rep <- apply(usageMat, 2, geoMean)
                names(rep) <- colnames(usageMat)

                idx <- .getIVYGAPidx(names(rep))
                ctIdx <- idx[["CT"]]
                itIdx <- idx[["IT"]]
                leIdx <- idx[["LE"]]
                mvpIdx <- idx[["MVP"]]
                panIdx <- idx[["PAN"]]
                
                ct <- rep[ctIdx]; ct <- ct[order(ct)]; ctF <- ecdf(rep[ctIdx])
                it <- rep[itIdx]; it <- it[order(it)]; itF <- ecdf(rep[itIdx])
                le <- rep[leIdx]; le <- le[order(le)]; leF <- ecdf(rep[leIdx])
                mvp <- rep[mvpIdx]; mvp <- mvp[order(mvp)]; mvpF <- ecdf(rep[mvpIdx])
                pan <- rep[panIdx]; pan <- pan[order(pan)]; panF <- ecdf(rep[panIdx])
                
                plot(ct, ctF(ct), col="red", xlab="GeoMean(GEPs)", ylab="ECDF", main=paste0("Community #", comm), type="l", lwd=3, xlim=c(0, max(rep)), ylim=c(0, 1))
                lines(it, itF(it), col="orange", lwd=3)
                lines(le, leF(le), col="gold", lwd=3)
                lines(mvp, mvpF(mvp), col="darkslategray4", lwd=3)
                lines(pan, panF(pan), col="burlywood4", lwd=3)
                
                legend("bottomright", 
                        legend=c(
                                paste0("CT (", length(ctIdx), ")"),  
                                paste0("IT (", length(ctIdx), ")"), 
                                paste0("LE (", length(ctIdx), ")"), 
                                paste0("MVP (", length(ctIdx), ")"), 
                                paste0("PAN (", length(ctIdx), ")")
                        ),
                        col = c("red", "orange", "gold", "darkslategray4", "burlywood4"), 
                        pch = 15, pt.cex = 2.8,
                        bty="n", ncol=1
                )
        } else {
                message(paste0("No metagenes in Community #", comm))
        }
}
dev.off()


# Get metagene-associated genes (MAGs)
expr <- readRDS(file.path(download_dir, "expr_GLASS.RDS"))
mags <- getMAGs(soObj, list(GLASS = expr))
message("MAGs for the metagenes: ", names(mags[["GLASS"]][[1]])[1])
mags[["GLASS"]][[1]][[1]] # [[data/cohort]][[community #]][[metagenes]]

# Get highly contributing genes
highlyContributingCommunityGenes <- contributingCommunityGenes(soObj, mags)
message("Highly contributing genes in Community #1 from GLASS:")
highlyContributingCommunityGenes[["GLASS"]][[1]]

q("no")
