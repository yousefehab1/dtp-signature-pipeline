test_that("attach_subtypes joins subtype calls and deduplicates properly", {
  clinical <- data.frame(
    Sample_ID = c("S1", "S2", "S3"),
    Age = c(50, 60, 70),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(
    Sample_ID = c("S1", "S2", "S2", "S3"),
    CMS = c("CMS1", "CMS2", "CMS2", "CMS3"),
    stringsAsFactors = FALSE
  )

  res <- attach_subtypes(clinical, subtypes, clinical_key = "Sample_ID")
  expect_equal(nrow(res), 3)
  expect_equal(res$CMS, c("CMS1", "CMS2", "CMS3"))

  # Null / empty subtype data returns original clinical unchanged
  expect_equal(attach_subtypes(clinical, NULL), clinical)
  expect_equal(attach_subtypes(clinical, data.frame()), clinical)
})

test_that("attach_subtypes supports custom key_fun for barcode-to-patient mapping", {
  clinical <- data.frame(
    Patient_ID = c("TCGA-AA-0001", "TCGA-AA-0002"),
    Stage = c("I", "II"),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(
    Sample_ID = c("TCGA-AA-0001-01A-01R", "TCGA-AA-0002-01A-01R"),
    PDS = c("PDS1", "PDS2"),
    stringsAsFactors = FALSE
  )

  res <- attach_subtypes(clinical, subtypes, clinical_key = "Patient_ID", key_fun = tcga_patient_id)
  expect_equal(nrow(res), 2)
  expect_equal(res$PDS, c("PDS1", "PDS2"))
  expect_false("Sample_ID" %in% colnames(res))
})

test_that("subtype_cohorts splits into named cohort data frames", {
  clinical <- data.frame(
    Sample_ID = paste0("S", 1:6),
    CMS = c("CMS1", "CMS2", "CMS1", "CMS3", NA, ""),
    stringsAsFactors = FALSE
  )

  cohorts <- subtype_cohorts(clinical, subtype_col = "CMS", prefix = "GSE")
  expect_setequal(names(cohorts), c("GSE_CMS1", "GSE_CMS2", "GSE_CMS3"))
  expect_equal(nrow(cohorts$GSE_CMS1), 2)
  expect_equal(nrow(cohorts$GSE_CMS2), 1)
  expect_equal(nrow(cohorts$GSE_CMS3), 1)

  # Non-existent subtype column returns empty list
  expect_equal(subtype_cohorts(clinical, "NonExistent", "GSE"), list())
})

test_that("build_symbol_matrix forces id_type='symbol' regardless of cfg", {
  testthat::skip_if_not_installed("hgu133plus2.db")
  testthat::skip_if_not_installed("SummarizedExperiment")

  cfg <- dtp_config(overrides = list(id_type = "ensembl"))

  # Microarray dispatch
  mat_micro <- matrix(c(5.0, 6.0), nrow = 1, dimnames = list("1007_s_at", c("S1", "S2")))
  res_micro <- build_symbol_matrix(cfg, mat_micro, source = "microarray")
  expect_equal(rownames(res_micro), "DDR1")

  # TCGA dispatch
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm_unstrand = matrix(c(10, 20, 30, 40), nrow = 2,
                                        dimnames = list(NULL, c("TCGA-AA-0001-01A", "TCGA-AA-0002-01A")))),
    rowData = S4Vectors::DataFrame(gene_id = c("ENSG00000001", "ENSG00000002"),
                                   gene_name = c("GENEA", "GENEB")),
    colData = S4Vectors::DataFrame(sample_type = c("Primary Tumor", "Primary Tumor"))
  )
  res_tcga <- build_symbol_matrix(cfg, se, source = "tcga")
  expect_equal(rownames(res_tcga), c("GENEA", "GENEB"))
})

test_that("call_cms and call_pds require realistic whole-transcriptome input", {
  # CMScaller requires matching hundreds of biological marker genes from its
  # classifier templates (or it halts with '<2 matched features/class').
  # PDSclassifier requires human gene sets for internal GSVA ssGSEA scoring.
  testthat::skip("call_cms/call_pds require full-transcriptome expression matrices with biological gene symbols")
})
