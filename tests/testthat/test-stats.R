test_that("metric_cols returns correct column mappings for OS and RFS", {
  expect_equal(metric_cols("OS"), c(t = "OS3Y_delay", e = "OS3Y_event"))
  expect_equal(metric_cols("RFS"), c(t = "RFS3Y_delay", e = "RFS3Y_event"))
})

test_that("get_wilcox_stats gates on min_group_n and non-numeric score", {
  cfg <- dtp_config(overrides = list(min_group_n = 5))

  # Non-numeric score column
  df_non_num <- data.frame(score = c("a", "b", "c", "d"), group = c("A", "A", "B", "B"))
  expect_equal(get_wilcox_stats(df_non_num, "score", "group", cfg),
               list(p = NA, r = NA, n = NA))

  # Group size smaller than min_group_n
  df_small <- data.frame(
    score = c(1, 2, 3, 4, 5, 6, 7),
    group = c("A", "A", "A", "A", "B", "B", "B") # A has 4, B has 3 < 5
  )
  expect_equal(get_wilcox_stats(df_small, "score", "group", cfg),
               list(p = NA, r = NA, n = NA))

  # Fewer than 2 groups
  df_one_grp <- data.frame(score = 1:10, group = rep("A", 10))
  expect_equal(get_wilcox_stats(df_one_grp, "score", "group", cfg),
               list(p = NA, r = NA, n = NA))
})

test_that("get_wilcox_stats calculates correct p-value and effect size", {
  cfg <- dtp_config(overrides = list(min_group_n = 3))
  g_low  <- c(1.2, 2.3, 3.1, 4.0)
  g_high <- c(5.5, 6.1, 7.8, 8.2)
  df <- data.frame(
    Score = c(g_low, g_high),
    Group = c(rep("Grp1", 4), rep("Grp2", 4)) # Grp1, Grp2 alphabetical
  )

  wt_ref <- suppressWarnings(stats::wilcox.test(g_low, g_high, exact = FALSE))
  res <- get_wilcox_stats(df, "Score", "Group", cfg)

  expect_equal(res$p, wt_ref$p.value)
  expect_equal(res$n, 8)
  expect_equal(res$r, as.numeric(1 - 2 * wt_ref$statistic / (4 * 4)))
})

test_that("get_kruskal_stats gates properly and computes epsilon-squared", {
  cfg <- dtp_config(overrides = list(min_group_n = 3))

  # 3 groups with >= 3 items each
  df <- data.frame(
    Score = c(1, 2, 3, 10, 11, 12, 20, 21, 22),
    Group = rep(c("A", "B", "C"), each = 3)
  )
  kt_ref <- stats::kruskal.test(Score ~ factor(Group), data = df)
  res <- get_kruskal_stats(df, "Score", "Group", cfg)

  expect_equal(res$p, kt_ref$p.value)
  expect_equal(res$n, 9)
  expect_equal(res$k, 3)
  n <- 9
  expect_equal(res$eps2, as.numeric(kt_ref$statistic) / ((n^2 - 1) / (n + 1)))

  # Filter drops groups with < min_group_n, leaving < 2 groups
  df_small_grps <- data.frame(
    Score = c(1, 2, 3, 10, 20),
    Group = c("A", "A", "A", "B", "C")
  )
  res_small <- get_kruskal_stats(df_small_grps, "Score", "Group", cfg)
  expect_true(is.na(res_small$p))
  expect_true(is.na(res_small$eps2))
})

test_that("get_km_pval gates on group size and event count", {
  cfg <- dtp_config(overrides = list(min_km_group_n = 3, min_events = 2))

  # Not enough events
  df_no_events <- data.frame(
    Score = 1:8,
    Time = c(10, 12, 14, 15, 20, 22, 25, 30),
    Event = c(1, 0, 0, 0, 0, 0, 0, 0) # Only 1 event < min_events (2)
  )
  expect_true(is.na(get_km_pval(df_no_events, "Score", "Time", "Event", cfg)))

  # Valid data with >= 2 events and >= 3 per KM group (median split)
  df_valid <- data.frame(
    Score = 1:8,
    Time = c(10, 12, 14, 15, 20, 22, 25, 30),
    Event = c(1, 1, 0, 0, 1, 1, 0, 0)
  )
  p <- get_km_pval(df_valid, "Score", "Time", "Event", cfg)
  expect_false(is.na(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("get_cox_stats gates on sample size and event count", {
  cfg <- dtp_config(overrides = list(min_cox_n = 6, min_events = 3))

  # Too few samples
  df_small <- data.frame(
    Score = 1:4,
    Time = c(10, 12, 14, 16),
    Event = c(1, 1, 1, 1)
  )
  expect_null(get_cox_stats(df_small, "Score", "Time", "Event", cfg))

  # Too few events
  df_few_events <- data.frame(
    Score = 1:8,
    Time = c(10, 12, 14, 16, 18, 20, 22, 24),
    Event = c(1, 1, 0, 0, 0, 0, 0, 0)
  )
  expect_null(get_cox_stats(df_few_events, "Score", "Time", "Event", cfg))

  # Valid data without complete separation
  df_valid <- data.frame(
    Score = c(0.1, 1.5, 0.4, 1.8, 0.2, 1.2, 0.8, 1.1, 0.5, 1.9),
    Time = c(5, 8, 10, 12, 15, 20, 25, 30, 35, 40),
    Event = c(1, 0, 1, 1, 0, 1, 0, 0, 1, 0)
  )
  res <- get_cox_stats(df_valid, "Score", "Time", "Event", cfg)
  expect_false(is.null(res))
  expect_true(is.numeric(res$HR))
  expect_true(is.numeric(res$P))
  expect_equal(res$N, 10)
  expect_equal(res$N_events, 5)
})

test_that("harmonize_crc_modifiers correctly unifies GSE39582 and TCGA-COAD", {
  # GSE39582
  gse_df <- data.frame(
    TNM_stage = c("1", "2", "3", "4", "NA"),
    MMR = c("pMMR", "dMMR", "pMMR", "dMMR", "pMMR"),
    CMS = c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"),
    PDS = c("PDS1", "PDS2", "PDS3", "PDS1", "Mixed"),
    stringsAsFactors = FALSE
  )
  gse_harm <- harmonize_crc_modifiers(gse_df, cohort = "GSE39582")
  expect_equal(levels(gse_harm$Stage_bin), c("Early", "Late"))
  expect_equal(as.character(gse_harm$Stage_bin), c("Early", "Early", "Late", "Late", NA))
  expect_equal(as.character(gse_harm$MSI_group), c("MSS", "MSI", "MSS", "MSI", "MSS"))
  expect_equal(levels(gse_harm$CMS)[1], "CMS2") # ref level
  expect_true(is.na(gse_harm$CMS[5]))

  # TCGA-COAD
  tcga_df <- data.frame(
    Stage = c("I", "II", "III", "IV", "X"),
    paper_MSI_status = c("MSS", "MSI-L", "MSI-H", "MSI-L", "MSS"),
    stringsAsFactors = FALSE
  )
  tcga_harm <- harmonize_crc_modifiers(tcga_df, cohort = "TCGA-COAD")
  expect_equal(as.character(tcga_harm$Stage_bin), c("Early", "Early", "Late", "Late", NA))
  expect_equal(as.character(tcga_harm$MSI_group), c("MSS", "MSS", "MSI", "MSS", "MSS"))
})

test_that("get_cox_adjusted and get_cox_interaction evaluate multivariable models", {
  set.seed(42)
  n <- 50
  cfg <- dtp_config(overrides = list(min_cox_n = 20, min_events = 10, min_epv = 2))
  df <- data.frame(
    Time = stats::rexp(n, rate = 0.05) + 1,
    Event = stats::rbinom(n, 1, 0.6),
    Score = stats::rnorm(n),
    Stage = factor(sample(c("Early", "Late"), n, replace = TRUE), levels = c("Early", "Late")),
    Age = stats::rnorm(n, mean = 60, sd = 10)
  )

  adj <- get_cox_adjusted(df, "Score", "Time", "Event", covars = c("Stage", "Age"), cfg = cfg)
  expect_false(is.null(adj))
  expect_true("HR" %in% names(adj))
  expect_true("LRT_P" %in% names(adj))
  expect_true("PH_P" %in% names(adj))
  expect_true("Delta_logHR" %in% names(adj))

  inter <- get_cox_interaction(df, "Score", "Time", "Event", modifier = "Stage", cfg = cfg)
  expect_false(is.null(inter))
  expect_true("Interaction_P" %in% names(inter))
  expect_equal(inter$K, 2)
  expect_true(is.data.frame(inter$per_level))
})

test_that("apply_fdr groups by cohort vs pancan configurations and errors on invalid family", {
  cfg <- dtp_config()

  # Construct a synthetic stats table with 2 datasets, 2 projects, 2 tests, 2 metrics
  df <- expand.grid(
    Dataset = c("D1", "D2"),
    Project = c("P1", "P2"),
    Test = c("Wilcoxon", "KM"),
    Metric = c("OS", "RFS"),
    stringsAsFactors = FALSE
  )
  df$Family <- c(rep("FamA", 8), rep("FamB", 8))
  df$Raw_P <- c(0.01, 0.04, 0.03, 0.02, 0.05, 0.001, 0.02, 0.04,
                0.01, 0.04, 0.03, 0.02, 0.05, 0.001, 0.02, 0.04)

  # cohort family groups by c("Dataset", "Project", "Test", "Metric") -> 1 row per group -> FDR = Raw_P
  res_cohort <- apply_fdr(df, cfg, family = "cohort")
  expect_true("FDR_P" %in% colnames(res_cohort))
  expect_true("Is_Significant" %in% colnames(res_cohort))
  expect_equal(res_cohort$FDR_P, df$Raw_P)

  # pancan family groups by c("Family", "Test") -> 4 rows per group -> FDR adjusted across 4 p-values
  res_pancan <- apply_fdr(df, cfg, family = "pancan")
  expect_false(isTRUE(all.equal(res_pancan$FDR_P, df$Raw_P)))

  # Invalid family argument errors out
  expect_error(apply_fdr(df, cfg, family = "invalid"), "should be one of")
})

test_that("run_survival_block runs end-to-end on synthetic cohort list", {
  cfg <- dtp_config(overrides = list(min_group_n = 3, min_km_group_n = 3, min_cox_n = 5, min_events = 2))

  c1 <- data.frame(
    OS3Y_delay = c(10, 20, 30, 36, 36, 15),
    OS3Y_event = c(1, 1, 0, 0, 0, 1),
    RFS3Y_delay = c(8, 15, 25, 36, 36, 12),
    RFS3Y_event = c(1, 1, 0, 0, 0, 1),
    Surv_3yr = factor(c("Dead_3yr", "Dead_3yr", "Alive_3yr", "Alive_3yr", "Alive_3yr", "Dead_3yr")),
    Recurrence_3yr = factor(c("Recurred", "Recurred", "Recurrence-Free", "Recurrence-Free", "Recurrence-Free", "Recurred")),
    SigA = c(1.2, 0.5, -0.3, -1.0, -0.8, 1.5),
    SigB = c(-0.5, 0.2, 1.1, 0.8, -0.2, 0.4)
  )

  cohorts <- list(CohortA = c1)
  res <- run_survival_block(cohorts, dataset = "TestDS", score_cols = c("SigA", "SigB"), cfg = cfg)

  # 1 cohort * 2 scores * 2 metrics (OS, RFS) * 3 tests (Wilcoxon, KM, Cox) = 12 rows
  expect_equal(nrow(res), 12)
  expect_setequal(res$Test, c("Wilcoxon", "KM", "Cox"))
  expect_setequal(res$Metric, c("OS", "RFS"))
  expect_setequal(res$Score, c("SigA", "SigB"))
})
