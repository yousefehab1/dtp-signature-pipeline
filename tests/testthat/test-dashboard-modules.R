test_that("dashboard module UI functions return valid shiny tags", {
  # Source dashboard modules if not in package namespace
  modules <- c(
    "mod_overview.R", "mod_crc.R", "mod_pancan.R",
    "mod_mets.R", "mod_tables.R", "mod_about.R"
  )
  for (m in modules) {
    p <- file.path("../../dashboard/R", m)
    if (file.exists(p)) source(p, local = FALSE)
    p2 <- file.path("dashboard/R", m)
    if (file.exists(p2)) source(p2, local = FALSE)
  }

  ui_funcs <- list(
    overview = mod_overview_ui,
    crc      = mod_crc_ui,
    pancan   = mod_pancan_ui,
    mets     = mod_mets_ui,
    tables   = mod_tables_ui,
    about    = mod_about_ui
  )

  for (nm in names(ui_funcs)) {
    fn <- ui_funcs[[nm]]
    ui_out <- fn(paste0("test_", nm))
    expect_true(inherits(ui_out, "shiny.tag") || inherits(ui_out, "shiny.tag.list"),
                info = paste("UI failed for module:", nm))
  }
})

test_that("dashboard module servers execute without error on populated bundle", {
  fake_bundle_pop <- list(
    gse = data.frame(
      Sample_ID = "S1", Up_ssGSEA = 0.5, Down_ssGSEA = -0.2,
      Surv_3yr = factor("Alive_3yr"), Recurrence_3yr = factor("Recurrence-Free"),
      OS3Y_delay = 36, OS3Y_event = 0, RFS3Y_delay = 36, RFS3Y_event = 0, Chemo_adj = "Y"
    ),
    tcga = data.frame(
      ID = "T1", Up_ssGSEA = 0.3, Down_ssGSEA = 0.1,
      Surv_3yr = factor("Alive_3yr"), Recurrence_3yr = factor("Recurrence-Free"),
      OS3Y_delay = 36, OS3Y_event = 0, RFS3Y_delay = 36, RFS3Y_event = 0
    ),
    panel_tbl = data.frame(
      Signature_ID = c("Up", "Down"),
      Signature_Name = c("DTP Up", "DTP Down"),
      Gene_ID = c("G1", "G2"),
      Source_Citation = c("Paper A", "Paper B"),
      stringsAsFactors = FALSE
    ),
    composite_defs = data.frame(
      Composite_ID = "Composite", Display_Name = "DTP Composite",
      Positive_Signature = "Up", Negative_Signature = "Down",
      Is_Default = TRUE, stringsAsFactors = FALSE
    ),
    stats_df = data.frame(
      Dataset = "GSE39582", Score = "Up_ssGSEA", Metric = "OS",
      Test = "Wilcoxon", Raw_P = 0.05, FDR_P = 0.05, stringsAsFactors = FALSE
    ),
    cfg = list(
      score_suffix = "_ssGSEA", os_cutpoint = 36, rfs_cutpoint = 36,
      min_group_n = 10, min_km_group_n = 10, min_events = 5
    ),
    built_at = Sys.time(),
    panel_mtime = Sys.time(),
    pancan_stats_df = data.frame(Project = "TCGA-COAD", Score = "Up_ssGSEA", HR_perSD = 1.2, Cox_FDR_P = 0.05),
    pancan_table_tbl = data.frame(Project = "TCGA-COAD", Score = "Up_ssGSEA", HR_perSD = 1.2, Cox_FDR_P = 0.05),
    treated_stats_df = data.frame(Project = "TCGA-COAD", Score = "Up_ssGSEA", HR_perSD = 1.4, Cox_FDR_P = 0.08),
    mets_gsea_pn = data.frame(ID = "Up", NES = 1.5, p.adjust = 0.01),
    mets_gsea_pm_adj = data.frame(ID = "Up", NES = 1.3, p.adjust = 0.02),
    mets_pvm_cmp = data.frame(ID = "Up", NES_unadjusted = 1.8, NES_adjusted = 1.3),
    mets_liver_df = data.frame(Sample = "S1", Score = 0.5),
    mets_roc_summary = data.frame(ID = "Up", AUC = 0.82),
    mets_ssgsea_pvalues = data.frame(Signature = "Up", p_adj_BH = 0.01),
    fig_dir = "output/figures",
    excel_path = "output/DTP_Results.xlsx"
  )

  # Test mod_overview_server
  expect_no_error({
    shiny::testServer(mod_overview_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })

  # Test mod_crc_server
  expect_no_error({
    shiny::testServer(mod_crc_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })

  # Test mod_pancan_server
  expect_no_error({
    shiny::testServer(mod_pancan_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })

  # Test mod_mets_server
  expect_no_error({
    shiny::testServer(mod_mets_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })

  # Test mod_tables_server
  expect_no_error({
    shiny::testServer(mod_tables_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })

  # Test mod_about_server
  expect_no_error({
    shiny::testServer(mod_about_server, args = list(bundle = fake_bundle_pop), {
      session$flushReact()
    })
  })
})

test_that("dashboard module servers execute without error on NULL/missing bundle", {
  # Test mod_overview_server with NULL
  expect_no_error({
    shiny::testServer(mod_overview_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })

  # Test mod_crc_server with NULL
  expect_no_error({
    shiny::testServer(mod_crc_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })

  # Test mod_pancan_server with NULL
  expect_no_error({
    shiny::testServer(mod_pancan_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })

  # Test mod_mets_server with NULL
  expect_no_error({
    shiny::testServer(mod_mets_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })

  # Test mod_tables_server with NULL
  expect_no_error({
    shiny::testServer(mod_tables_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })

  # Test mod_about_server with NULL
  expect_no_error({
    shiny::testServer(mod_about_server, args = list(bundle = NULL), {
      session$flushReact()
    })
  })
})
