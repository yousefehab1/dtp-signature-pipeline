test_that("canonical panel loads with the expected signatures", {
  cfg <- dtp_config(root = test_path("..", ".."))
  panel <- load_signature_panel(cfg)

  expect_setequal(unique(panel$Signature_ID),
                   c("Up", "Down", "Fetal", "revSC", "RSC", "CBC", "MYC", "IBD"))
  expect_false(anyDuplicated(panel[c("Signature_ID", "Gene_ID")]) > 0)
})

test_that("panel_to_list drops NA/blank genes and matches panel row counts", {
  cfg <- dtp_config(root = test_path("..", ".."))
  panel <- load_signature_panel(cfg)
  panel_list <- panel_to_list(panel, id_type = "ensembl")

  expect_equal(sum(lengths(panel_list)), nrow(panel))
  expect_true(all(vapply(panel_list, function(g) all(!is.na(g) & trimws(g) != ""), logical(1))))
})

test_that("load_signature_panel rejects a panel with duplicate rows", {
  cfg <- dtp_config(root = test_path("..", ".."))
  bad <- tempfile(fileext = ".csv")
  on.exit(unlink(bad))
  writeLines(c("Signature_ID,Gene_ID", "Up,ENSG1", "Up,ENSG1"), bad)

  expect_error(load_signature_panel(cfg, path = bad), "duplicate")
})

test_that("check_updown_overlap warns when Up and Down share a gene", {
  panel_list <- list(Up = c("A", "B"), Down = c("B", "C"))
  expect_warning(check_updown_overlap(panel_list, "Up", "Down"), "1 gene")
})

test_that("panel_hash is stable under row reordering", {
  panel <- data.frame(Signature_ID = c("Up", "Down"), Gene_ID = c("G1", "G2"))
  panel_reordered <- panel[c(2, 1), ]
  expect_equal(panel_hash(panel), panel_hash(panel_reordered))
})

test_that("load_composite_defs exposes the default Up-Down composite", {
  cfg <- dtp_config(root = test_path("..", ".."))
  cdefs <- load_composite_defs(cfg)

  expect_true(any(cdefs$Is_Default))
  default_row <- cdefs[cdefs$Is_Default, ]
  expect_equal(default_row$Positive_Signature[1], "Up")
  expect_equal(default_row$Negative_Signature[1], "Down")
})
