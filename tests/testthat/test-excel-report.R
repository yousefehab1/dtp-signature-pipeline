test_that("build_excel_workbook creates valid workbook with CRC sheets and Session_Info in order", {
  cfg <- dtp_config()
  tmp_xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp_xlsx), add = TRUE)

  panel_tbl <- data.frame(
    Signature_ID = c("Up", "Down"),
    Signature_Name = c("DTP Up", "DTP Down"),
    Gene_ID = c("G1", "G2"),
    stringsAsFactors = FALSE
  )
  composite_defs <- data.frame(
    Composite_ID = "Composite",
    Display_Name = "DTP Composite",
    Positive_Signature = "Up",
    Negative_Signature = "Down",
    Is_Default = TRUE,
    stringsAsFactors = FALSE
  )
  crc <- list(
    stats_df = data.frame(Dataset = "GSE39582", FDR_P = c(0.01, 0.2)),
    subtype_survival_tbl = data.frame(Subgroup = "CMS1", HR_perSD = 1.5, FDR_P = 0.04),
    adj_df = data.frame(Model = "Clinicopath", FDR_P = 0.02),
    int_df = data.frame(Modifier = "CMS", FDR_P = 0.5),
    subtype_stats = data.frame(Subtype_Axis = "CMS", FDR_P = 0.001)
  )

  out <- build_excel_workbook(cfg, tmp_xlsx, panel_tbl, composite_defs, crc)

  expect_true(file.exists(tmp_xlsx))
  sheets <- openxlsx::getSheetNames(tmp_xlsx)
  expected_sheets <- c(
    "README",
    "Signature_Panel",
    "Composite_Definitions",
    "CRC_Univariable_Survival",
    "CRC_Subgroup_Survival",
    "CRC_Adjusted_Cox",
    "CRC_Interaction_Cox",
    "CRC_Subtype_Characterization",
    "Session_Info"
  )
  expect_equal(sheets, expected_sheets)

  # Check Session_Info sheet content
  sess_df <- openxlsx::read.xlsx(tmp_xlsx, sheet = "Session_Info")
  expect_gt(nrow(sess_df), 10)
  expect_true(any(grepl("R Version", sess_df$Property)))
  expect_true(any(grepl("dtpsig", sess_df$Property)))
  expect_true(any(grepl("survival", sess_df$Property)))
})

test_that("build_excel_workbook works with all-empty results lists", {
  cfg <- dtp_config()
  tmp_xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp_xlsx), add = TRUE)

  panel_tbl <- data.frame(
    Signature_ID = "Up",
    Signature_Name = "DTP Up",
    Gene_ID = "G1",
    stringsAsFactors = FALSE
  )
  composite_defs <- data.frame(
    Composite_ID = "Composite",
    Display_Name = "DTP Composite",
    Positive_Signature = "Up",
    Negative_Signature = "None",
    Is_Default = TRUE,
    stringsAsFactors = FALSE
  )

  out <- build_excel_workbook(cfg, tmp_xlsx, panel_tbl, composite_defs, crc = list(), pancan = list(), mets = list())
  expect_true(file.exists(tmp_xlsx))
  sheets <- openxlsx::getSheetNames(tmp_xlsx)
  expect_true("Session_Info" %in% sheets)
  expect_true("README" %in% sheets)
  expect_true("Signature_Panel" %in% sheets)
  expect_true("Composite_Definitions" %in% sheets)
})
