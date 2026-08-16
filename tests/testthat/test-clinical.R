test_that("derive_endpoints derives landmark endpoints for all censoring and event cases", {
  # Synthetic data covering the 4 cases + NA
  # 1. Event before cutpoint (20 mo < 36 mo)
  # 2. Event after cutpoint (50 mo > 36 mo)
  # 3. Censored before cutpoint (20 mo < 36 mo)
  # 4. Censored after cutpoint (50 mo > 36 mo)
  # 5. Missing values
  df <- data.frame(
    OS_event = c(1, 1, 0, 0, NA),
    OS_delay = c(20, 50, 20, 50, NA),
    Recurrence_event = c(1, 1, 0, 0, NA),
    Recurrence_delay = c(20, 50, 20, 50, NA)
  )

  res <- derive_endpoints(df,
                          os_event = "OS_event", os_time = "OS_delay",
                          rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
                          os_cut = 36, rfs_cut = 36)

  # Check Surv_3yr levels and values
  expect_equal(levels(res$Surv_3yr), c("Dead_3yr", "Alive_3yr"))
  expect_equal(as.character(res$Surv_3yr), c("Dead_3yr", "Alive_3yr", "Alive_3yr", "Alive_3yr", NA))

  # Check Recurrence_3yr levels and values
  expect_equal(levels(res$Recurrence_3yr), c("Recurred", "Recurrence-Free"))
  expect_equal(as.character(res$Recurrence_3yr), c("Recurred", "Recurrence-Free", NA, "Recurrence-Free", NA))

  # Check OS 3Y numeric columns
  expect_equal(res$OS3Y_event, c(1, 0, 0, 0, NA))
  expect_equal(res$OS3Y_delay, c(20, 36, 20, 36, NA))

  # Check RFS 3Y numeric columns
  expect_equal(res$RFS3Y_event, c(1, 0, 0, 0, NA))
  expect_equal(res$RFS3Y_delay, c(20, 36, 20, 36, NA))
})

test_that("load_tcga_cdr parses CDR CSV correctly", {
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv))

  cdr_raw <- data.frame(
    bcr_patient_barcode = c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"),
    type = c("COAD", "COAD", "READ"),
    OS = c(1, 0, 1),
    OS.time = c(304.4, 608.8, 913.2),
    PFI = c(1, 0, 0),
    PFI.time = c(152.2, 608.8, 913.2),
    ajcc_pathologic_tumor_stage = c("Stage IIA", "Stage III", "Stage I"),
    stringsAsFactors = FALSE
  )
  write.csv(cdr_raw, tmp_csv, row.names = FALSE)

  cfg <- dtp_config(overrides = list(cdr_file = tmp_csv, cdr_rfs_type = "PFI"))
  cdr <- load_tcga_cdr(cfg)

  expect_equal(nrow(cdr), 3)
  expect_setequal(colnames(cdr), c("bcr_patient_barcode", "type", "OS", "OS.time",
                                   "PFI", "PFI.time", "ajcc_pathologic_tumor_stage",
                                   "Patient_ID", "Project_ID", "Stage",
                                   "OS_event", "OS_months", "DFS_event", "DFS_months"))
  expect_equal(cdr$Patient_ID, c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"))
  expect_equal(cdr$Project_ID, c("TCGA-COAD", "TCGA-COAD", "TCGA-READ"))
  expect_equal(cdr$Stage, c("II", "III", "I"))
  expect_equal(cdr$OS_months, c(10, 20, 30))
  expect_equal(cdr$DFS_months, c(5, 20, 30))
})
