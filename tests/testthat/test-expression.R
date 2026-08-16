test_that("gene_set_coverage calculates percentage of genes present in matrix", {
  mat <- matrix(1:12, nrow = 4, dimnames = list(c("G1", "G2", "G3", "G4"), c("S1", "S2", "S3")))

  expect_equal(gene_set_coverage(mat, c("G1", "G2", "G5", "G6")), 50)
  expect_equal(gene_set_coverage(mat, c("G1", "G2", "G3", "G4")), 100)
  expect_equal(gene_set_coverage(mat, c("G7", "G8")), 0)
})

test_that("prep_microarray_symbols maps probes and collapses duplicates", {
  testthat::skip_if_not_installed("hgu133plus2.db")
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("limma")

  cfg_sym <- dtp_config(overrides = list(id_type = "symbol"))
  # 1007_s_at is DDR1 on HG-U133 Plus 2.0
  mat <- matrix(c(5.0, 7.0, 6.0, 8.0), nrow = 2, byrow = TRUE,
                dimnames = list(c("1007_s_at", "1007_s_at"), c("S1", "S2")))

  res <- prep_microarray_symbols(cfg_sym, mat)
  expect_equal(rownames(res), "DDR1")
  expect_equal(res["DDR1", "S1"], 5.5) # Average of 5.0 and 6.0
  expect_equal(res["DDR1", "S2"], 7.5) # Average of 7.0 and 8.0

  cfg_ens <- dtp_config(overrides = list(id_type = "ensembl"))
  res_ens <- prep_microarray_symbols(cfg_ens, mat)
  expect_true(grepl("^ENSG", rownames(res_ens)))
})

test_that("prep_tcga_tpm log2-transforms, strips Ensembl versions, and filters sample types", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("limma")

  cfg_ens <- dtp_config(overrides = list(
    id_type = "ensembl",
    tumour_codes = c("Primary Tumor")
  ))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm_unstrand = matrix(
      c(3, 7,
        15, 31,
        1, 1),
      nrow = 2, ncol = 3,
      dimnames = list(NULL, c("S1_Tumor", "S2_Tumor", "S3_Normal"))
    )),
    rowData = S4Vectors::DataFrame(
      gene_id = c("ENSG00000141510.16", "ENSG00000133703.12"),
      gene_name = c("TP53", "KRAS")
    ),
    colData = S4Vectors::DataFrame(
      sample_type = c("Primary Tumor", "Primary Tumor", "Solid Tissue Normal")
    )
  )

  res_ens <- prep_tcga_tpm(cfg_ens, se)
  # Filtered out normal sample (S3) -> only 2 samples
  expect_equal(colnames(res_ens), c("S1_Tumor", "S2_Tumor"))
  expect_equal(rownames(res_ens), c("ENSG00000141510", "ENSG00000133703"))
  # S1_Tumor: TP53=3 -> log2(3+1)=2, KRAS=7 -> log2(7+1)=3
  # S2_Tumor: TP53=15 -> log2(15+1)=4, KRAS=31 -> log2(31+1)=5
  expect_equal(res_ens["ENSG00000141510", "S1_Tumor"], 2)
  expect_equal(res_ens["ENSG00000133703", "S1_Tumor"], 3)
  expect_equal(res_ens["ENSG00000141510", "S2_Tumor"], 4)
  expect_equal(res_ens["ENSG00000133703", "S2_Tumor"], 5)

  # Symbol target mode
  cfg_sym <- dtp_config(overrides = list(
    id_type = "symbol",
    tumour_codes = c("Primary Tumor")
  ))
  res_sym <- prep_tcga_tpm(cfg_sym, se)
  expect_equal(rownames(res_sym), c("TP53", "KRAS"))

  # Error on missing assay
  se_bad_assay <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, 2, 2))
  )
  expect_error(prep_tcga_tpm(cfg_ens, se_bad_assay), "Assay 'tpm_unstrand' not found")
})
