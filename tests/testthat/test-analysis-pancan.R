test_that(".pancan_is_transient_gdc_error accurately classifies errors", {
  # Transient errors
  e_timeout <- simpleError("Connection timed out after 10000 milliseconds")
  e_503     <- simpleError("HTTP 503 Service Unavailable")
  e_curl    <- simpleError("curl error 56: SSL read: error:00000000:lib(0):func(0):reason(0), errno 0")
  e_net     <- simpleError("Network connection was lost")

  expect_true(.pancan_is_transient_gdc_error(e_timeout))
  expect_true(.pancan_is_transient_gdc_error(e_503))
  expect_true(.pancan_is_transient_gdc_error(e_curl))
  expect_true(.pancan_is_transient_gdc_error(e_net))

  # Non-transient errors
  e_subscript <- structure(simpleError("Subscript out of bounds"), class = c("vctrs_error_subscript_oob", "error", "condition"))
  e_parse     <- simpleError("Parse error in column specification")
  e_type      <- simpleError("Cannot coerce type closure to character")

  expect_false(.pancan_is_transient_gdc_error(e_subscript))
  expect_false(.pancan_is_transient_gdc_error(e_parse))
  expect_false(.pancan_is_transient_gdc_error(e_type))
})

test_that(".pancan_batch_correct handles dynamic score columns and adjusts batches", {
  set.seed(42)
  n_per_proj <- 25
  # Synthetic master with 2 projects and 4 arbitrary score columns
  proj_A <- data.frame(
    Patient_ID = paste0("A_", seq_len(n_per_proj)),
    Project_ID = "TCGA-PA",
    Sig1_ssGSEA = rnorm(n_per_proj, mean = 2, sd = 0.5),
    Sig2_ssGSEA = rnorm(n_per_proj, mean = 5, sd = 0.5),
    Sig3_ssGSEA = rnorm(n_per_proj, mean = -1, sd = 0.5),
    Sig4_ssGSEA = rnorm(n_per_proj, mean = 0.5, sd = 0.5),
    stringsAsFactors = FALSE
  )
  proj_B <- data.frame(
    Patient_ID = paste0("B_", seq_len(n_per_proj)),
    Project_ID = "TCGA-PB",
    Sig1_ssGSEA = rnorm(n_per_proj, mean = 8, sd = 0.5), # intentional batch shift
    Sig2_ssGSEA = rnorm(n_per_proj, mean = 1, sd = 0.5),
    Sig3_ssGSEA = rnorm(n_per_proj, mean = 3, sd = 0.5),
    Sig4_ssGSEA = rnorm(n_per_proj, mean = -2, sd = 0.5),
    stringsAsFactors = FALSE
  )
  master <- rbind(proj_A, proj_B)
  score_cols <- c("Sig1_ssGSEA", "Sig2_ssGSEA", "Sig3_ssGSEA", "Sig4_ssGSEA")

  bc_res <- .pancan_batch_correct(master, score_cols)

  expect_equal(bc_res$corr_cols, paste0(score_cols, "_Corrected"))
  expect_equal(nrow(bc_res$master), 50)
  expect_true(all(bc_res$corr_cols %in% colnames(bc_res$master)))

  # Before batch correction, means across projects differed significantly
  mean_A_raw <- mean(master$Sig1_ssGSEA[master$Project_ID == "TCGA-PA"])
  mean_B_raw <- mean(master$Sig1_ssGSEA[master$Project_ID == "TCGA-PB"])
  expect_true(abs(mean_A_raw - mean_B_raw) > 5)

  # After batch correction, means across projects are aligned
  mean_A_corr <- mean(bc_res$master$Sig1_ssGSEA_Corrected[bc_res$master$Project_ID == "TCGA-PA"])
  mean_B_corr <- mean(bc_res$master$Sig1_ssGSEA_Corrected[bc_res$master$Project_ID == "TCGA-PB"])
  expect_true(abs(mean_A_corr - mean_B_corr) < 1e-10)
})

test_that(".pancan_stats_block runs statistics across dynamic scores and applies pancan FDR", {
  set.seed(42)
  cfg <- dtp_config(overrides = list(min_group_n = 5, min_cox_n = 10, min_events = 3))

  # Synthetic master table with survival endpoints
  n <- 60
  master <- data.frame(
    Patient_ID       = paste0("P", seq_len(n)),
    Project_ID       = rep(c("TCGA-P1", "TCGA-P2"), each = n / 2),
    SigA_ssGSEA      = rnorm(n),
    SigB_ssGSEA      = rnorm(n),
    OS_event         = rbinom(n, 1, 0.4),
    OS_delay         = runif(n, 5, 50),
    Recurrence_event = rbinom(n, 1, 0.3),
    Recurrence_delay = runif(n, 5, 50),
    stringsAsFactors = FALSE
  )
  master <- derive_endpoints(master, os_event = "OS_event", os_time = "OS_delay",
                             rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
                             os_cut = 36, rfs_cut = 36)

  score_cols <- c("SigA_ssGSEA", "SigB_ssGSEA")
  bc <- .pancan_batch_correct(master, score_cols)

  stats_res <- .pancan_stats_block(bc$master, score_cols, bc$corr_cols, cfg)
  stats_df  <- stats_res$stats_df
  cox_sum   <- stats_res$cox_summary

  expect_true(is.data.frame(stats_df))
  expect_true(all(c("Dataset", "Project", "Test", "Metric", "Score", "Raw_P", "Family", "FDR_P") %in% colnames(stats_df)))

  # Check family assignment
  expect_true(all(stats_df$Family[stats_df$Project == "Pan-Cancer"] == "PanCancer"))
  expect_true(all(stats_df$Family[stats_df$Project != "Pan-Cancer"] == "PerCohort"))

  # Check dynamic scores are evaluated
  expect_true(all(c("SigA_ssGSEA_Corrected", "SigB_ssGSEA_Corrected") %in% stats_df$Score[stats_df$Project == "Pan-Cancer"]))
  expect_true(all(c("SigA_ssGSEA", "SigB_ssGSEA") %in% stats_df$Score[stats_df$Project != "Pan-Cancer"]))

  # Check Cox summary
  expect_true(all(!is.na(cox_sum$HR)))
  expect_true(all(c("HR_lower", "HR_upper", "FDR_P") %in% colnames(cox_sum)))
})

test_that(".pancan_treated_stats_block gates cohorts below min_cox_n for aggregate", {
  set.seed(42)
  cfg <- dtp_config(overrides = list(min_cox_n = 15, min_events = 3, min_group_n = 5))

  # Project 1 has 20 treated, Project 2 has 8 treated (< min_cox_n 15)
  n1 <- 20
  n2 <- 8
  p1 <- data.frame(
    Patient_ID = paste0("T1_", seq_len(n1)), Project_ID = "TCGA-LARGE",
    Score1_ssGSEA = rnorm(n1), OS_event = rbinom(n1, 1, 0.4), OS_delay = runif(n1, 5, 50),
    Recurrence_event = rbinom(n1, 1, 0.3), Recurrence_delay = runif(n1, 5, 50),
    stringsAsFactors = FALSE
  )
  p2 <- data.frame(
    Patient_ID = paste0("T2_", seq_len(n2)), Project_ID = "TCGA-SMALL",
    Score1_ssGSEA = rnorm(n2), OS_event = rbinom(n2, 1, 0.4), OS_delay = runif(n2, 5, 50),
    Recurrence_event = rbinom(n2, 1, 0.3), Recurrence_delay = runif(n2, 5, 50),
    stringsAsFactors = FALSE
  )
  mt_all <- rbind(p1, p2)
  mt_all <- derive_endpoints(mt_all, os_event = "OS_event", os_time = "OS_delay",
                             rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
                             os_cut = 36, rfs_cut = 36)

  counts <- data.frame(
    Project_ID  = c("TCGA-LARGE", "TCGA-SMALL"),
    N_total     = c(25, 10),
    N_treated   = c(20, 8),
    N_untreated = c(3, 1),
    N_unknown   = c(2, 1),
    stringsAsFactors = FALSE
  )

  score_cols <- "Score1_ssGSEA"
  # Since only 1 cohort meets >= 15 treated, aggregate should be skipped (< 2 cohorts)
  stats_df <- .pancan_treated_stats_block(mt_all, counts, score_cols, cfg)

  expect_true(is.data.frame(stats_df))
  expect_false("Pan-Cancer" %in% stats_df$Project)
  expect_true("TCGA-LARGE" %in% stats_df$Project)
  expect_true("TCGA-SMALL" %in% stats_df$Project)
})

test_that(".pancan_treatment_flags caching semantics: loads from cache and does not cache errors", {
  tmp_cache <- file.path(tempdir(), "test_treat_cache")
  dir.create(tmp_cache, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  cfg <- dtp_config(overrides = list(cache_dir = tmp_cache))

  # 1. Direct cache hit from a pre-populated cache_dir
  fake_df <- data.frame(
    Patient_ID = c("TCGA-01-0001", "TCGA-01-0002"),
    Project_ID = "TCGA-TEST1",
    pharma_yes = c(TRUE, FALSE), pharma_no = c(FALSE, TRUE),
    rad_yes = c(FALSE, FALSE), rad_no = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  saveRDS(fake_df, file.path(tmp_cache, "treat_flags_TCGA_TEST1.rds"))

  res_hit <- .pancan_treatment_flags(cfg, "TCGA-TEST1")
  expect_equal(res_hit, fake_df)

  # 2. Failed network pull is NOT cached
  testthat::with_mocked_bindings(
    GDCquery_clinic = function(project, type) stop("Simulated network timeout"),
    .package = "TCGAbiolinks",
    {
      res_fail <- .pancan_treatment_flags(cfg, "TCGA-FAILPROJ")
      expect_equal(nrow(res_fail), 0)
      expect_false(file.exists(file.path(tmp_cache, "treat_flags_TCGA_FAILPROJ.rds")))
    }
  )

  # 3. Successful network pull IS cached
  fake_gdc_clinic <- data.frame(
    submitter_id = c("TCGA-02-0001", "TCGA-02-0002"),
    treatments_pharmaceutical_treatment_or_therapy = c("yes", "no"),
    treatments_radiation_treatment_or_therapy = c("no", "yes"),
    stringsAsFactors = FALSE
  )
  testthat::with_mocked_bindings(
    GDCquery_clinic = function(project, type) fake_gdc_clinic,
    .package = "TCGAbiolinks",
    {
      res_succ <- .pancan_treatment_flags(cfg, "TCGA-SUCCPROJ")
      expect_equal(nrow(res_succ), 2)
      expect_equal(res_succ$Patient_ID, c("TCGA-02-0001", "TCGA-02-0002"))
      expect_true(res_succ$pharma_yes[1])
      expect_true(res_succ$rad_yes[2])
      expect_true(file.exists(file.path(tmp_cache, "treat_flags_TCGA_SUCCPROJ.rds")))
    }
  )
})
