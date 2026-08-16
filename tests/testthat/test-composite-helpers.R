test_that(".crc_modules dynamically derives labels and order from panel and composite defs", {
  panel_tbl <- data.frame(
    Signature_ID = c("SigA", "SigA", "SigB"),
    Signature_Name = c("Signature Alpha", "Signature Alpha", "Signature Beta"),
    Gene_ID = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  composite_defs <- data.frame(
    Composite_ID = c("Comp1", "Comp2"),
    Display_Name = c("Composite One", "Composite Two"),
    Positive_Signature = c("SigA", "SigB"),
    Negative_Signature = c("SigB", "SigA"),
    Is_Default = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  cfg <- list(score_suffix = "_ssGSEA")

  # Test standard set including default composite, panel sigs, non-default composite, and novel score
  score_cols <- c("SigB_ssGSEA", "Novel_score", "Comp1_ssGSEA", "SigA_ssGSEA", "Comp2_ssGSEA")
  res <- .crc_modules(score_cols, panel_tbl, composite_defs, cfg)

  # Default composite Comp1 comes first, then panel sigs in panel order (SigA, SigB), then Comp2, then Novel_score
  expect_equal(res$cols, c("Comp1_ssGSEA", "SigA_ssGSEA", "SigB_ssGSEA", "Comp2_ssGSEA", "Novel_score"))
  expect_equal(res$labels[["Comp1_ssGSEA"]], "Composite One")
  expect_equal(res$labels[["SigA_ssGSEA"]], "Signature Alpha")
  expect_equal(res$labels[["SigB_ssGSEA"]], "Signature Beta")
  expect_equal(res$labels[["Comp2_ssGSEA"]], "Composite Two")
  expect_equal(res$labels[["Novel_score"]], "Novel_score") # Not dropped, not error
  expect_equal(res$levels, c("Composite One", "Signature Alpha", "Signature Beta", "Composite Two", "Novel_score"))
})

test_that(".subgroup_category and .subgroup_label parse known project formats", {
  # Category classification
  expect_equal(.subgroup_category("GSE_CMS1"), "CMS")
  expect_equal(.subgroup_category("TCGA_CMS4"), "CMS")
  expect_equal(.subgroup_category("GSE_PDS2"), "PDS")
  expect_equal(.subgroup_category("TCGA_Mixed"), "PDS")
  expect_equal(.subgroup_category("GSE_Stage1"), "Stage")
  expect_equal(.subgroup_category("TCGA_Stage4"), "Stage")
  expect_equal(.subgroup_category("GSE_All_MSS"), "MSI")
  expect_equal(.subgroup_category("TCGA_All_MSI"), "MSI")
  expect_true(is.na(.subgroup_category("GSE_Stage3_Treated")))
  expect_true(is.na(.subgroup_category("GSE_All_Patients")))

  # Label formatting
  expect_equal(.subgroup_label("GSE_CMS1"), "CMS1")
  expect_equal(.subgroup_label("TCGA_Stage1"), "Stage I")
  expect_equal(.subgroup_label("GSE_Stage3"), "Stage III")
  expect_equal(.subgroup_label("TCGA_All_MSS"), "All MSS")
})

test_that(".subtype_pair_stats calculates pairwise Mann-Whitney and effect size", {
  # Synthetic clinical dataframe with 2 CMS levels
  set.seed(42)
  clin <- data.frame(
    CMS = c(rep("CMS1", 15), rep("CMS2", 15), rep("CMS3", 15), rep("CMS4", 15)),
    Up_ssGSEA = c(rnorm(15, mean = 1.0), rnorm(15, mean = 0.5), rnorm(15, mean = 0.8), rnorm(15, mean = 1.2)),
    stringsAsFactors = FALSE
  )
  cfg <- list(min_group_n = 5)

  res <- .subtype_pair_stats(clin, dataset = "GSE39582", axis = "CMS", score = "Up_ssGSEA", cfg = cfg)

  expect_equal(nrow(res), 6) # 4 choose 2 = 6 pairs
  expect_true(all(c("Dataset", "Subtype_Axis", "Score", "Group_A", "Group_B", "Effect_r", "Raw_P", "Is_Testable") %in% colnames(res)))
  expect_true(all(res$Is_Testable))
  expect_false(any(is.na(res$Raw_P)))
  expect_false(any(is.na(res$Effect_r)))

  # Check signed effect r: CMS1 vs CMS2 has higher mean in CMS1 -> negative r
  r_12 <- res$Effect_r[res$Group_A == "CMS1" & res$Group_B == "CMS2"]
  expect_true(r_12 < 0)
})
