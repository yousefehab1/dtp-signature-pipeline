test_that("detect_id_type identifies ensembl vs symbol correctly", {
  expect_equal(detect_id_type(c("ENSG00000141510", "ENSG00000133703", "ENSG00000123456")), "ensembl")
  expect_equal(detect_id_type(c("TP53", "KRAS", "EGFR", "CDX2")), "symbol")
  expect_equal(detect_id_type(character(0)), "symbol")
  expect_equal(detect_id_type(c(NA_character_, "")), "symbol")
  # 80% threshold test: 4 out of 5 ENSG (80%) -> not > 0.8 -> symbol; 5 out of 5 -> ensembl
  expect_equal(detect_id_type(c("ENSG1", "ENSG2", "ENSG3", "ENSG4", "TP53")), "symbol")
  expect_equal(detect_id_type(c("ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5", "TP53")), "ensembl")
})

test_that("strip_ensembl_version removes .version suffix", {
  ids <- c("ENSG00000141510.16", "ENSG00000133703.12", "ENSG00000123456", "TP53.1")
  expect_equal(strip_ensembl_version(ids),
               c("ENSG00000141510", "ENSG00000133703", "ENSG00000123456", "TP53"))
})

test_that("convert_gene_ids converts symbols to ensembl via org.Hs.eg.db", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("AnnotationDbi")

  cfg <- dtp_config()
  converted <- convert_gene_ids(cfg, c("TP53", "KRAS"), from = "SYMBOL", to = "ENSEMBL")
  expect_true("ENSG00000141510" %in% converted) # TP53
  expect_true("ENSG00000133703" %in% converted) # KRAS
})

test_that("harmonise_matrix_ids harmonises row names to target id_type", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("AnnotationDbi")

  # Convert symbol matrix to ensembl target
  cfg_ens <- dtp_config(overrides = list(id_type = "ensembl"))
  mat_sym <- matrix(c(10, 20, 30, 40), nrow = 2,
                    dimnames = list(c("TP53", "KRAS"), c("S1", "S2")))
  mat_harm_ens <- harmonise_matrix_ids(cfg_ens, mat_sym, current_type = "symbol")
  expect_true("ENSG00000141510" %in% rownames(mat_harm_ens))
  expect_true("ENSG00000133703" %in% rownames(mat_harm_ens))

  # Target ensembl when already ensembl with version numbers
  mat_ver <- matrix(c(10, 20), nrow = 1,
                    dimnames = list("ENSG00000141510.16", c("S1", "S2")))
  mat_harm_ver <- harmonise_matrix_ids(cfg_ens, mat_ver, current_type = "ensembl")
  expect_equal(rownames(mat_harm_ver), "ENSG00000141510")

  # Convert ensembl matrix to symbol target
  cfg_sym <- dtp_config(overrides = list(id_type = "symbol"))
  mat_ens <- matrix(c(10, 20), nrow = 1,
                    dimnames = list("ENSG00000141510", c("S1", "S2")))
  mat_harm_sym <- harmonise_matrix_ids(cfg_sym, mat_ens, current_type = "ensembl")
  expect_equal(rownames(mat_harm_sym), "TP53")
})
