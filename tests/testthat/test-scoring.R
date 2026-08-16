test_that("add_composite computes positive minus negative score with correct column name", {
  df <- data.frame(
    Up_ssGSEA = c(0.8, -0.2, 0.5),
    Down_ssGSEA = c(-0.1, 0.4, 0.1)
  )

  res <- add_composite(df, positive = "Up", negative = "Down", name = "Composite", suffix = "_ssGSEA")
  expect_equal(res$Composite_ssGSEA, c(0.9, -0.6, 0.4))

  # If one of the columns is missing, data frame remains unchanged
  res_missing <- add_composite(df, positive = "NonExistent", negative = "Down", name = "Comp")
  expect_false("Comp_ssGSEA" %in% colnames(res_missing))
})

test_that("add_all_composites loops over definitions table", {
  df <- data.frame(
    Up_ssGSEA = c(1.0, 2.0),
    Down_ssGSEA = c(0.5, 1.0),
    SigA_ssGSEA = c(0.2, 0.4),
    SigB_ssGSEA = c(0.1, 0.2)
  )

  cdefs <- data.frame(
    Composite_ID = c("Comp1", "Comp2"),
    Display_Name = c("Composite 1", "Composite 2"),
    Positive_Signature = c("Up", "SigA"),
    Negative_Signature = c("Down", "SigB"),
    Is_Default = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  res <- add_all_composites(df, cdefs, suffix = "_ssGSEA")
  expect_equal(res$Comp1_ssGSEA, c(0.5, 1.0))
  expect_equal(res$Comp2_ssGSEA, c(0.1, 0.2))
})

test_that("tcga_patient_id extracts 12-char patient barcode accurately", {
  barcodes <- c(
    "TCGA-AA-3873-01A-01R-0905-07",
    "tcga-a6-2675-01a",
    "TCGA.AZ.6598.01A",
    "GSM123456",
    "",
    NA
  )

  expected <- c(
    "TCGA-AA-3873",
    "TCGA-A6-2675",
    "TCGA-AZ-6598",
    NA_character_,
    NA_character_,
    NA_character_
  )

  expect_equal(tcga_patient_id(barcodes), expected)
})

test_that("run_ssgsea scores expression matrix with GeneSetCollection", {
  testthat::skip_if_not_installed("GSVA")
  testthat::skip_if_not_installed("GSEABase")

  mat <- matrix(
    c(10, 8, 2, 1,
      9, 7, 3, 2,
      1, 2, 8, 9),
    nrow = 4, byrow = FALSE,
    dimnames = list(c("G1", "G2", "G3", "G4"), c("S1", "S2", "S3"))
  )

  gs <- GSEABase::GeneSetCollection(list(
    GSEABase::GeneSet(c("G1", "G2"), setName = "SigUp"),
    GSEABase::GeneSet(c("G3", "G4"), setName = "SigDown")
  ))

  res <- run_ssgsea(mat, gs, suffix = "_ssGSEA")
  expect_equal(nrow(res), 3)
  expect_setequal(colnames(res), c("SigUp_ssGSEA", "SigDown_ssGSEA"))
  expect_true(is.numeric(res$SigUp_ssGSEA))
  expect_true(is.numeric(res$SigDown_ssGSEA))
})
