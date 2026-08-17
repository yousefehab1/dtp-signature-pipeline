test_that(".pancan_rescale correctly calculates per-1-SD HR and handles non-convergent fits", {
  # 1. Toy HR + SD hand-computed check
  d_toy <- data.frame(
    HR = 2.0,
    HR_lower = 1.5,
    HR_upper = 2.5,
    stringsAsFactors = FALSE
  )
  sd_val <- 0.5
  res_toy <- .pancan_rescale(d_toy, sd_val)

  expect_equal(res_toy$hr_sd, sqrt(2.0))
  expect_equal(res_toy$lo_sd, sqrt(1.5))
  expect_equal(res_toy$hi_sd, sqrt(2.5))
  expect_equal(res_toy$dir, "Higher hazard (HR>1)")
  expect_false(res_toy$offscale)
  expect_false(res_toy$nonconv)
  expect_equal(res_toy$hr_d, sqrt(2.0))

  # 2. Non-convergent fits (Inf / 0 / NaN HR) must not crash and must not clamp to boundary
  d_nonconv <- data.frame(
    HR = c(Inf, 0, 1.5),
    HR_lower = c(0.1, 0, 1.1),
    HR_upper = c(Inf, 1.2, 2.0),
    stringsAsFactors = FALSE
  )
  res_nc <- .pancan_rescale(d_nonconv, c(1.0, 1.0, 1.0))

  # First two rows are non-convergent
  expect_true(res_nc$nonconv[1])
  expect_true(res_nc$nonconv[2])
  expect_false(res_nc$nonconv[3])

  # Display coords for non-convergent rows must be NA_real_ (not squished to 10 or 0.1)
  expect_true(is.na(res_nc$hr_d[1]))
  expect_true(is.na(res_nc$hr_d[2]))
  expect_equal(res_nc$hr_d[3], 1.5)

  expect_true(is.na(res_nc$dir[1]))
  expect_true(is.na(res_nc$dir[2]))
  expect_equal(res_nc$dir[3], "Higher hazard (HR>1)")
})

test_that(".pancan_cohort_order and build_pancan_forest produce expected cohort row order", {
  # Synthetic master table with 3 cohorts
  master <- data.frame(
    Project_ID = rep(c("TCGA-A", "TCGA-B", "TCGA-C"), each = 20),
    Up_ssGSEA = rep(c(1.0, 1.0, 1.0), each = 20) + rnorm(60, sd = 0.5),
    OS3Y_event = rep(1, 60),
    OS3Y_delay = rep(20, 60),
    RFS3Y_event = rep(1, 60),
    RFS3Y_delay = rep(20, 60),
    stringsAsFactors = FALSE
  )

  # Stats table where TCGA-C has lowest HR, TCGA-A middle, TCGA-B highest
  stats_df <- data.frame(
    Family = "PerCohort",
    Test = "Cox",
    Score = "Up_ssGSEA",
    Metric = rep(c("OS", "RFS"), 3),
    Project = rep(c("TCGA-A", "TCGA-B", "TCGA-C"), each = 2),
    N = 20,
    N_events = 10,
    HR = c(1.2, 1.1, 2.5, 2.0, 0.6, 0.7),
    HR_lower = c(1.0, 0.9, 1.8, 1.5, 0.4, 0.5),
    HR_upper = c(1.5, 1.3, 3.5, 2.8, 0.9, 0.95),
    Raw_P = 0.01,
    FDR_P = 0.02,
    Is_Testable = TRUE,
    Is_Significant = TRUE,
    stringsAsFactors = FALSE
  )

  ord <- .pancan_cohort_order(stats_df, master, order_score = "Up_ssGSEA")
  # Expected ascending order by OS per-SD HR: TCGA-C, TCGA-A, TCGA-B
  expect_equal(ord, c("TCGA-C", "TCGA-A", "TCGA-B"))

  # build_pancan_forest returns a ggplot object with these factor levels
  p <- build_pancan_forest(stats_df, master, score = "Up_ssGSEA", score_lab = "DTP Up", order_score = "Up_ssGSEA")
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$Project), c("TCGA-C", "TCGA-A", "TCGA-B"))
})

test_that("build_pancan_aggregate, build_pancan_table, and build_pancan_violins generate correct outputs", {
  set.seed(42)
  n <- 60
  master <- data.frame(
    Patient_ID = paste0("P", seq_len(n)),
    Project_ID = rep(c("TCGA-A", "TCGA-B"), each = n / 2),
    Up_ssGSEA = rnorm(n, mean = 1),
    Down_ssGSEA = rnorm(n, mean = -1),
    Composite_ssGSEA = rnorm(n, mean = 0.5),
    Up_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    Down_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    Composite_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    OS3Y_event = rbinom(n, 1, 0.4),
    OS3Y_delay = runif(n, 5, 50),
    RFS3Y_event = rbinom(n, 1, 0.3),
    RFS3Y_delay = runif(n, 5, 50),
    Surv_3yr = sample(c("Dead_3yr", "Alive_3yr"), n, replace = TRUE),
    Recurrence_3yr = sample(c("Recurred", "Recurrence-Free"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  stats_df <- data.frame(
    Family = c(rep("PerCohort", 12), rep("PanCancer", 6), rep("PanCancer", 6)),
    Test = c(rep("Cox", 12), rep("Cox", 6), rep("Wilcoxon", 6)),
    Score = c(
      rep(c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA"), each = 4),
      rep(c("Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"), each = 2),
      rep(c("Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"), each = 2)
    ),
    Metric = c(
      rep(rep(c("OS", "RFS"), each = 2), 3),
      rep(c("OS", "RFS"), 3),
      rep(c("OS", "RFS"), 3)
    ),
    Project = c(
      rep(c("TCGA-A", "TCGA-B"), 6),
      rep("Pan-Cancer", 6),
      rep("Pan-Cancer", 6)
    ),
    N = 30,
    N_events = 10,
    HR = 1.3,
    HR_lower = 1.05,
    HR_upper = 1.6,
    C_index = 0.55,
    Effect_r = -0.15,
    Raw_P = 0.02,
    FDR_P = 0.03,
    Is_Testable = TRUE,
    Is_Significant = TRUE,
    stringsAsFactors = FALSE
  )

  cfg <- dtp_config(root = test_path("..", ".."))
  composite_defs <- load_composite_defs(cfg)

  # 1. Aggregate forest
  p_agg <- build_pancan_aggregate(stats_df, master, cfg = cfg, composite_defs = composite_defs)
  expect_s3_class(p_agg, "ggplot")

  # 2. Pancan table (Table 3.5)
  tbl <- build_pancan_table(stats_df, master, cfg = cfg, composite_defs = composite_defs)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(colnames(tbl), c(
    "Project", "Score", "Endpoint", "N", "N_events", "Is_Testable",
    "HR_perSD", "HR_lower", "HR_upper", "C_index", "Cox_FDR_P", "Is_Significant"
  ))
  expect_true(any(tbl$Project == "Pan-Cancer (batch-corrected)"))

  # 3. Pancan violins (Fig 7)
  res_violins <- build_pancan_violins(stats_df, master, cfg = cfg, composite_defs = composite_defs)
  expect_type(res_violins, "list")
  expect_s3_class(res_violins$plot, "ggplot")
  expect_type(res_violins$caption, "character")
  expect_true(grepl("Figure 7", res_violins$caption))

  # 4. Fig 6 composite
  pA <- build_pancan_forest(stats_df, master, score = "Up_ssGSEA", score_lab = "DTP Up")
  pB <- build_pancan_forest(stats_df, master, score = "Composite_ssGSEA", score_lab = "DTP Composite")
  pC <- build_pancan_forest(stats_df, master, score = "Down_ssGSEA", score_lab = "DTP Down")
  p_comp <- build_fig6_pancan_survival_composite(pA, pB, pC)
  expect_s3_class(p_comp, "patchwork")
})

test_that("build_pancan_composites handles execution, saves files, and isolates group failures via tryCatch", {
  set.seed(42)
  n <- 60
  master <- data.frame(
    Patient_ID = paste0("P", seq_len(n)),
    Project_ID = rep(c("TCGA-COAD", "TCGA-BRCA"), each = n / 2),
    Up_ssGSEA = rnorm(n, mean = 1),
    Down_ssGSEA = rnorm(n, mean = -1),
    Composite_ssGSEA = rnorm(n, mean = 0.5),
    Up_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    Down_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    Composite_ssGSEA_Corrected = rnorm(n, mean = 0, sd = 0.5),
    OS3Y_event = rbinom(n, 1, 0.4),
    OS3Y_delay = runif(n, 5, 50),
    RFS3Y_event = rbinom(n, 1, 0.3),
    RFS3Y_delay = runif(n, 5, 50),
    Surv_3yr = sample(c("Dead_3yr", "Alive_3yr"), n, replace = TRUE),
    Recurrence_3yr = sample(c("Recurred", "Recurrence-Free"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  stats_df <- data.frame(
    Family = c(rep("PerCohort", 12), rep("PanCancer", 6), rep("PanCancer", 6)),
    Test = c(rep("Cox", 12), rep("Cox", 6), rep("Wilcoxon", 6)),
    Score = c(
      rep(c("Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA"), each = 4),
      rep(c("Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"), each = 2),
      rep(c("Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"), each = 2)
    ),
    Metric = c(
      rep(rep(c("OS", "RFS"), each = 2), 3),
      rep(c("OS", "RFS"), 3),
      rep(c("OS", "RFS"), 3)
    ),
    Project = c(
      rep(c("TCGA-COAD", "TCGA-BRCA"), 6),
      rep("Pan-Cancer", 6),
      rep("Pan-Cancer", 6)
    ),
    N = 30,
    N_events = 10,
    HR = 1.3,
    HR_lower = 1.05,
    HR_upper = 1.6,
    C_index = 0.55,
    Effect_r = -0.15,
    Raw_P = 0.02,
    FDR_P = 0.03,
    Is_Testable = TRUE,
    Is_Significant = TRUE,
    stringsAsFactors = FALSE
  )

  pancan_result <- list(
    stats_df = stats_df,
    master = master
  )

  cfg <- dtp_config(root = test_path("..", ".."))
  tmp_dir <- tempfile("pancan_fig_test_")
  dir.create(tmp_dir, recursive = TRUE)

  # Full run should succeed on all groups
  res <- build_pancan_composites(pancan_result, cfg = cfg, out_dir = tmp_dir)
  expect_type(res, "list")
  expect_true(all(unlist(res$status)))
  expect_s3_class(res$pancan_table_tbl, "tbl_df")
  expect_true(file.exists(file.path(tmp_dir, "Fig6A_pancan_forest_Up.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig6B_pancan_forest_Composite.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig6C_pancan_forest_Down.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig6_pancan_survival_composite.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig_Pancan_aggregate_forest.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig7_pancan_violins.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "Fig6_pancan_survival_caption.txt")))
  expect_true(file.exists(file.path(tmp_dir, "Fig7_pancan_violins_caption.txt")))

  # Broken input isolation test: corrupt Surv_3yr so violins fails while other groups succeed
  broken_master <- master
  broken_master$Surv_3yr <- NULL
  broken_res_input <- list(stats_df = stats_df, master = broken_master)

  tmp_dir_broken <- tempfile("pancan_broken_test_")
  dir.create(tmp_dir_broken, recursive = TRUE)

  suppressWarnings({
    res_broken <- build_pancan_composites(broken_res_input, cfg = cfg, out_dir = tmp_dir_broken)
  })

  # Violins fails, but forest plots (A, B, C, composite), aggregate, and table still succeed
  expect_false(res_broken$status[["violins (Fig7)"]])
  expect_true(res_broken$status[["A (Up forest)"]])
  expect_true(res_broken$status[["B (Composite forest)"]])
  expect_true(res_broken$status[["C (Down forest)"]])
  expect_true(res_broken$status[["composite"]])
  expect_true(res_broken$status[["aggregate (appendix)"]])
  expect_true(res_broken$status[["table"]])
})

test_that("build_excel_workbook creates workbook with PanCancer sheets when pancan result passed", {
  cfg <- dtp_config(root = test_path("..", ".."))
  tmp_xlsx <- tempfile(fileext = ".xlsx")

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
  pancan_data <- list(
    pancan_table_tbl = data.frame(
      Project = c("TCGA-COAD", "Pan-Cancer (batch-corrected)"),
      Score = c("Up_ssGSEA", "Up_ssGSEA_Corrected"),
      Endpoint = c("OS", "OS"),
      HR_perSD = c(1.35, 1.09),
      FDR_P = c(0.04, 0.001),
      stringsAsFactors = FALSE
    ),
    treated_stats_df = data.frame(
      Project = "TCGA-COAD",
      Score = "Up_ssGSEA",
      Metric = "OS",
      HR = 1.25,
      FDR_P = 0.08,
      stringsAsFactors = FALSE
    )
  )

  build_excel_workbook(cfg, tmp_xlsx, panel_tbl, composite_defs, crc = list(), pancan = pancan_data)

  expect_true(file.exists(tmp_xlsx))
  sheets <- openxlsx::getSheetNames(tmp_xlsx)
  expect_true("PanCancer_Survival" %in% sheets)
  expect_true("PanCancer_Treated" %in% sheets)
  expect_true("Session_Info" %in% sheets)
  expect_equal(length(sheets), 11)
})
