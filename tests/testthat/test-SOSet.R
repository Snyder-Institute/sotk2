test_that("SOSet rejects invalid NMFobjL / NMFrankL / corMet inputs", {
        pvObjFile <- system.file("extdata", "PV.RDS", package = "sotk2")
        skip_if(pvObjFile == "", "PV.RDS not installed.")

        pvNmfObj <- readRDS(pvObjFile)
        dataL <- list(PV = pvNmfObj)
        rankL <- list(PV = c(2:6))

        # NMFobjL passed as the raw NMF object (not a list) — names mismatch
        expect_error(
                SOSet(NMFobjL = pvNmfObj, NMFrankL = rankL,
                      dataCol = c("PV" = "#7570B3"), corMet = "spearman"),
                "Check the names in NMFrankL"
        )

        # Unsupported correlation method
        expect_error(
                SOSet(NMFobjL = dataL, NMFrankL = rankL,
                      dataCol = c("PV" = "#7570B3"), corMet = "concordance_index"),
                "Provide either pearson, spearman, or kendall"
        )
})

test_that("SOSet constructs successfully on the bundled PV demo", {
        pvObjFile <- system.file("extdata", "PV.RDS", package = "sotk2")
        skip_if(pvObjFile == "", "PV.RDS not installed.")

        pvNmfObj <- readRDS(pvObjFile)
        soSet <- SOSet(
                NMFobjL  = list(PV = pvNmfObj),
                NMFrankL = list(PV = c(2:6)),
                dataCol  = c("PV" = "#7570B3"),
                corMet   = "spearman"
        )
        expect_s4_class(soSet, "SOSet")
        expect_true(isSymmetric(soSet@corMat))
})
