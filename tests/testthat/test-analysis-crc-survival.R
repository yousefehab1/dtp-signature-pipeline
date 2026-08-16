test_that(".crc_check_gse39582_consistency detects inconsistent label prefixes", {
  # Consistent characteristics columns
  clinical_good <- data.frame(
    characteristics_ch1.2 = c("Sex: Male", "Sex: Female"),
    characteristics_ch1.3 = c("age.at.diagnosis (year): 50", "age.at.diagnosis (year): 60"),
    stringsAsFactors = FALSE
  )
  expect_true(.crc_check_gse39582_consistency(clinical_good))

  # Inconsistent label in characteristics_ch1.2
  clinical_bad <- data.frame(
    characteristics_ch1.2 = c("Sex: Male", "age.at.diagnosis (year): 50"),
    characteristics_ch1.3 = c("age.at.diagnosis (year): 50", "age.at.diagnosis (year): 60"),
    stringsAsFactors = FALSE
  )
  expect_error(
    .crc_check_gse39582_consistency(clinical_bad),
    "characteristics_ch1.2 has 2 labels"
  )
})

test_that(".crc_parse_numeric_field parses numbers and drops unparseable values with message", {
  # All valid numbers
  res_good <- .crc_parse_numeric_field(c("age: 45", "age: 62.5", "age: 10"), "Age")
  expect_equal(res_good, c(45, 62.5, 10))

  # Contains unparseable strings
  expect_message(
    res_mixed <- .crc_parse_numeric_field(c("age: 45", "age: unknown", "age: N/A"), "Age"),
    "Age: 2 value\\(s\\) unparseable"
  )
  expect_equal(res_mixed, c(45, NA_real_, NA_real_))
})

test_that(".crc_build_gse_cohorts partitions clinical data accurately", {
  clinical <- data.frame(
    Sample_ID = paste0("S", 1:5),
    Chemo_adj = c("Y", "Y", "N", "N", "Y"),
    TNM_stage = c("3", "4", "2", "1", "3"),
    MMR       = c("pMMR", "pMMR", "pMMR", "dMMR", "dMMR"),
    CMS       = c("CMS1", "CMS2", "CMS3", "CMS4", "CMS1"),
    PDS       = c("PDS1", "PDS2", "PDS3", "PDS1", "PDS2"),
    stringsAsFactors = FALSE
  )

  cohorts <- .crc_build_gse_cohorts(clinical)

  expect_equal(nrow(cohorts$GSE_All_Patients), 5)
  expect_equal(cohorts$GSE_Treated$Sample_ID, c("S1", "S2", "S5"))
  expect_equal(cohorts$GSE_Stage3_4_Treated$Sample_ID, c("S1", "S2", "S5"))
  expect_equal(cohorts$GSE_Stage3_Treated$Sample_ID, c("S1", "S5"))
  expect_equal(cohorts$GSE_Stage3_4_Treated_MSS$Sample_ID, c("S1", "S2"))

  expect_equal(cohorts$GSE_All_Treated$Sample_ID, c("S1", "S2", "S5"))
  expect_equal(cohorts$GSE_All_MSS$Sample_ID, c("S1", "S2", "S3"))
  expect_equal(cohorts$GSE_Stage2_Untreated_MSS$Sample_ID, "S3")
  expect_equal(cohorts$GSE_Stage3_Treated_MSS$Sample_ID, "S1")
  expect_equal(cohorts$GSE_Stage34_Treated_MSS$Sample_ID, c("S1", "S2"))

  expect_equal(cohorts$GSE_Stage1$Sample_ID, "S4")
  expect_equal(cohorts$GSE_Stage2$Sample_ID, "S3")
  expect_equal(cohorts$GSE_Stage3$Sample_ID, c("S1", "S5"))
  expect_equal(cohorts$GSE_Stage4$Sample_ID, "S2")
  expect_equal(cohorts$GSE_All_MSI$Sample_ID, c("S4", "S5"))

  expect_equal(cohorts$GSE_CMS1$Sample_ID, c("S1", "S5"))
  expect_equal(cohorts$GSE_CMS2$Sample_ID, "S2")
  expect_equal(cohorts$GSE_PDS1$Sample_ID, c("S1", "S4"))
  expect_equal(cohorts$GSE_PDS2$Sample_ID, c("S2", "S5"))
})

test_that(".crc_build_tcga_cohorts partitions TCGA clinical data accurately", {
  clinical_tcga <- data.frame(
    ID               = paste0("T", 1:4),
    Treatment_Status = c("Treated", "Not Treated", "Not Treated", "Treated"),
    Stage            = c("III", "I", "II", "IV"),
    paper_MSI_status = c("MSS", "MSI-H", "MSI-L", "MSS"),
    CMS              = c("CMS2", "CMS1", "CMS3", "CMS4"),
    PDS              = c("PDS1", "PDS2", "PDS3", "PDS1"),
    stringsAsFactors = FALSE
  )

  cohorts <- .crc_build_tcga_cohorts(clinical_tcga)

  expect_equal(nrow(cohorts$TCGA_All_Patients), 4)
  expect_equal(cohorts$TCGA_Untreated$ID, c("T2", "T3"))
  expect_equal(cohorts$TCGA_Treated$ID, c("T1", "T4"))
  expect_equal(cohorts$TCGA_Stage1_Untreated$ID, "T2")
  expect_equal(cohorts$TCGA_Stage2_Untreated$ID, "T3")
  expect_equal(cohorts$TCGA_Stage1$ID, "T2")
  expect_equal(cohorts$TCGA_Stage2$ID, "T3")
  expect_equal(cohorts$TCGA_Stage3$ID, "T1")
  expect_equal(cohorts$TCGA_Stage4$ID, "T4")
  expect_equal(cohorts$TCGA_All_MSS$ID, c("T1", "T4"))
  expect_equal(cohorts$TCGA_All_MSI$ID, c("T2", "T3"))
  expect_equal(cohorts$TCGA_CMS1$ID, "T2")
  expect_equal(cohorts$TCGA_PDS1$ID, c("T1", "T4"))
})

test_that("core_scores and modifier intersections degrade gracefully with missing entries", {
  cfg <- dtp_config()
  score_cols <- c("Up_ssGSEA", "Other_ssGSEA")
  core_scores <- intersect(paste0(c("Up", "Down", "Composite"), cfg$score_suffix), score_cols)
  expect_equal(core_scores, "Up_ssGSEA")

  # When no core scores match
  core_empty <- intersect(paste0(c("Up", "Down", "Composite"), cfg$score_suffix), c("A_ssGSEA", "B_ssGSEA"))
  expect_equal(core_empty, character(0))

  # Modifier intersection
  cd_cols <- c("CMS", "Age", "Stage_bin")
  mods <- intersect(cfg$crc_modifiers, cd_cols)
  expect_setequal(mods, c("CMS", "Stage_bin"))
})
