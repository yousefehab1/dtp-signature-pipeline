# ==============================================================================
# dashboard/R/mod_signature_explorer.R - Signature Explorer Shiny Module
#
# Interactive exploration of canonical and live custom signature pairs across
# CRC cohorts (GSE39582, TCGA-COAD) and clinical endpoints (OS, RFS).
# ==============================================================================

#' Signature Explorer UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_signature_explorer_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    # Staleness warning banner
    shiny::uiOutput(ns("staleness_alert")),

    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = "Query Parameters",
        width = 340,
        shiny::selectInput(
          ns("cohort"),
          label = "Cohort",
          choices = c(
            "GSE39582 (whole cohort)"   = "GSE_all",
            "GSE39582 (adjuvant-treated)" = "GSE_treated",
            "TCGA-COAD"                 = "TCGA_COAD"
          ),
          selected = "GSE_all"
        ),
        shiny::selectInput(
          ns("endpoint"),
          label = "Survival Endpoint",
          choices = c(
            "Overall Survival (OS)"        = "OS",
            "Recurrence-Free Survival (RFS)" = "RFS"
          ),
          selected = "OS"
        ),
        shiny::radioButtons(
          ns("score_mode"),
          label = "Score Mode",
          choices = c(
            "Predefined signature / composite" = "predefined",
            "Custom signature pair (live)"      = "custom"
          ),
          selected = "predefined"
        ),
        shiny::conditionalPanel(
          condition = "input.score_mode == 'predefined'",
          ns = ns,
          shiny::selectInput(
            ns("predefined_score"),
            label = "Predefined Signature / Composite",
            choices = c("Loading..." = "")
          )
        ),
        shiny::conditionalPanel(
          condition = "input.score_mode == 'custom'",
          ns = ns,
          shiny::selectInput(
            ns("custom_pos"),
            label = "Positive Signature (+)",
            choices = c("Loading..." = "")
          ),
          shiny::selectInput(
            ns("custom_neg"),
            label = "Negative Signature (-)",
            choices = c("None" = "None")
          ),
          shiny::helpText("Live score = Positive - Negative (or Positive if Negative is None). Recomputed on-the-fly.")
        )
      ),

      # Main Content Area
      bslib::layout_column_wrap(
        width = 1/3,
        bslib::value_box(
          title = "Wilcoxon Rank-Sum",
          value = shiny::textOutput(ns("vb_wilcox_val")),
          showcase = NULL,
          shiny::textOutput(ns("vb_wilcox_sub")),
          theme = "primary"
        ),
        bslib::value_box(
          title = "Kaplan-Meier (Log-Rank)",
          value = shiny::textOutput(ns("vb_km_val")),
          showcase = NULL,
          shiny::textOutput(ns("vb_km_sub")),
          theme = "primary"
        ),
        bslib::value_box(
          title = "Univariate Cox PH (per 1 SD)",
          value = shiny::textOutput(ns("vb_cox_val")),
          showcase = NULL,
          shiny::textOutput(ns("vb_cox_sub")),
          theme = "primary"
        )
      ),

      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header(shiny::strong("3-Year Outcome Distribution")),
          shiny::plotOutput(ns("plot_violin"), height = "420px")
        ),
        bslib::card(
          bslib::card_header(shiny::strong("Kaplan-Meier Survival Curve")),
          shiny::plotOutput(ns("plot_km"), height = "420px")
        )
      ),

      bslib::card(
        bslib::card_header(shiny::strong("Live Query Statistics Table")),
        DT::dataTableOutput(ns("stats_table")),
        bslib::card_footer(
          class = "text-muted small",
          "Live query - not FDR-corrected; see the static CRC report for BH-adjusted results across the full signature panel."
        )
      )
    )
  )
}

#' Signature Explorer Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @param is_stale Logical flag indicating if `panel.csv` is newer than bundle.
#' @export
mod_signature_explorer_server <- function(id, bundle, is_stale = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 1. Staleness banner
    output$staleness_alert <- shiny::renderUI({
      if (is.null(bundle)) {
        bslib::card(
          class = "bg-danger text-white mb-3",
          bslib::card_body(
            "⚠️ Results bundle missing! Please run `pipeline/build_results_bundle.R` to generate precomputed data."
          )
        )
      } else if (isTRUE(is_stale)) {
        bslib::card(
          class = "bg-warning text-dark mb-3",
          bslib::card_body(
            "⚠️ Notice: `inst/signatures/panel.csv` has changed since this bundle was built. Re-run `pipeline/build_results_bundle.R` to refresh precomputed results."
          )
        )
      } else {
        NULL
      }
    })

    # 2. Populate signature choices on load
    shiny::observe({
      shiny::req(bundle)
      suffix <- if (!is.null(bundle$cfg$score_suffix)) bundle$cfg$score_suffix else "_ssGSEA"
      score_cols <- unique(c(
        grep(paste0(suffix, "$"), names(bundle$gse), value = TRUE),
        grep(paste0(suffix, "$"), names(bundle$tcga), value = TRUE)
      ))

      mod_info <- dtpsig:::.crc_modules(score_cols, bundle$panel_tbl, bundle$composite_defs, bundle$cfg)
      pred_choices <- stats::setNames(mod_info$cols, unname(mod_info$labels[mod_info$cols]))

      def_sel <- if ("Composite_ssGSEA" %in% mod_info$cols) "Composite_ssGSEA" else mod_info$cols[1]
      shiny::updateSelectInput(session, "predefined_score", choices = pred_choices, selected = def_sel)

      sig_df <- unique(bundle$panel_tbl[, c("Signature_ID", "Signature_Name")])
      pos_choices <- stats::setNames(sig_df$Signature_ID, paste0(sig_df$Signature_Name, " (", sig_df$Signature_ID, ")"))
      neg_choices <- c("None" = "None", pos_choices)

      def_pos <- if ("Up" %in% sig_df$Signature_ID) "Up" else sig_df$Signature_ID[1]
      shiny::updateSelectInput(session, "custom_pos", choices = pos_choices, selected = def_pos)
      shiny::updateSelectInput(session, "custom_neg", choices = neg_choices, selected = "None")
    })

    # 3. Reactive computation pipeline
    reactive_res <- shiny::reactive({
      shiny::req(bundle)
      cohort_id <- input$cohort
      endpoint  <- input$endpoint
      mode      <- input$score_mode
      shiny::req(cohort_id, endpoint, mode)

      # A. Select clinical cohort
      if (cohort_id == "GSE_all") {
        clin <- bundle$gse
        cohort_display <- "GSE39582 (whole cohort)"
        ds_name <- "GSE39582"
      } else if (cohort_id == "GSE_treated") {
        clin <- bundle$gse[which(bundle$gse$Chemo_adj == "Y"), , drop = FALSE]
        cohort_display <- "GSE39582 (adjuvant-treated)"
        ds_name <- "GSE39582"
      } else {
        clin <- bundle$tcga
        cohort_display <- "TCGA-COAD"
        ds_name <- "TCGA-COAD"
      }

      # B. Compute or extract score
      suffix <- if (!is.null(bundle$cfg$score_suffix)) bundle$cfg$score_suffix else "_ssGSEA"
      score_cols <- unique(c(
        grep(paste0(suffix, "$"), names(bundle$gse), value = TRUE),
        grep(paste0(suffix, "$"), names(bundle$tcga), value = TRUE)
      ))
      mod_info <- dtpsig:::.crc_modules(score_cols, bundle$panel_tbl, bundle$composite_defs, bundle$cfg)

      if (mode == "predefined") {
        sc <- input$predefined_score
        shiny::req(sc)
        if (!sc %in% colnames(clin)) return(NULL)
        score_vec <- suppressWarnings(as.numeric(clin[[sc]]))
        score_label <- if (sc %in% names(mod_info$labels)) mod_info$labels[[sc]] else sc
      } else {
        pos_id <- input$custom_pos
        neg_id <- input$custom_neg
        shiny::req(pos_id)
        score_vec <- tryCatch(
          dtpsig:::.compute_live_score(clin, pos_id, neg_id, suffix = suffix),
          error = function(e) NULL
        )
        if (is.null(score_vec)) return(NULL)

        pos_name <- bundle$panel_tbl$Signature_Name[match(pos_id, bundle$panel_tbl$Signature_ID)]
        if (is.null(neg_id) || is.na(neg_id) || neg_id == "None" || neg_id == "") {
          score_label <- paste0(pos_name, " (", pos_id, ")")
        } else {
          neg_name <- bundle$panel_tbl$Signature_Name[match(neg_id, bundle$panel_tbl$Signature_ID)]
          score_label <- paste0(pos_id, " - ", neg_id, " (", pos_name, " - ", neg_name, ")")
        }
      }

      clin$live_score <- score_vec
      mc <- metric_cols(endpoint)
      group_col <- if (endpoint == "OS") "Surv_3yr" else "Recurrence_3yr"
      cfg <- bundle$cfg

      # C. Wilcoxon Rank-Sum test
      d_w <- clin[!is.na(clin[[group_col]]) & clin[[group_col]] != "" & !is.na(clin$live_score), , drop = FALSE]
      tbl_w <- table(droplevels(factor(d_w[[group_col]])))
      min_grp <- if (!is.null(cfg$min_group_n)) cfg$min_group_n else 10
      w_testable <- length(tbl_w) >= 2 && all(tbl_w >= min_grp)
      ws <- if (w_testable) get_wilcox_stats(clin, "live_score", group_col, cfg) else list(p = NA_real_, r = NA_real_, n = nrow(d_w))

      # D. Kaplan-Meier Log-Rank test
      d_km <- clin[!is.na(clin[[mc[["t"]]]]) & !is.na(clin[[mc[["e"]]]]) & !is.na(clin$live_score), , drop = FALSE]
      km_p <- get_km_pval(clin, "live_score", mc[["t"]], mc[["e"]], cfg)
      km_testable <- !is.na(km_p)

      # E. Univariate Cox PH test
      cs <- get_cox_stats(clin, "live_score", mc[["t"]], mc[["e"]], cfg)
      cox_testable <- !is.null(cs)

      sd_val <- if (nrow(d_km) >= 2) dtpsig:::.score_sd(d_km$live_score) else NA_real_
      hr_sd   <- if (cox_testable && is.finite(cs$HR) && !is.na(sd_val)) exp(log(cs$HR) * sd_val) else NA_real_
      hr_l_sd <- if (cox_testable && is.finite(cs$HR_lower) && !is.na(sd_val)) exp(log(cs$HR_lower) * sd_val) else NA_real_
      hr_u_sd <- if (cox_testable && is.finite(cs$HR_upper) && !is.na(sd_val)) exp(log(cs$HR_upper) * sd_val) else NA_real_

      list(
        clin           = clin,
        d_w            = d_w,
        d_km           = d_km,
        cohort_display = cohort_display,
        ds_name        = ds_name,
        endpoint       = endpoint,
        group_col      = group_col,
        t_col          = mc[["t"]],
        e_col          = mc[["e"]],
        score_label    = score_label,
        cfg            = cfg,
        ws             = ws,
        w_testable     = w_testable,
        km_p           = km_p,
        km_testable    = km_testable,
        cs             = cs,
        cox_testable   = cox_testable,
        hr_sd          = hr_sd,
        hr_l_sd        = hr_l_sd,
        hr_u_sd        = hr_u_sd,
        sd_val         = sd_val
      )
    })

    # 4. Value Boxes
    output$vb_wilcox_val <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res) || !res$w_testable || is.na(res$ws$p)) {
        "Gated / Untestable"
      } else {
        sprintf("r = %.2f (p = %.3g)", res$ws$r, res$ws$p)
      }
    })
    output$vb_wilcox_sub <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res)) return("")
      sprintf("n = %d across 3-yr outcome groups", res$ws$n)
    })

    output$vb_km_val <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res) || !res$km_testable || is.na(res$km_p)) {
        "Gated / Untestable"
      } else {
        sprintf("Log-rank p = %.3g", res$km_p)
      }
    })
    output$vb_km_sub <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res)) return("")
      n_ev <- if (nrow(res$d_km) > 0) sum(res$d_km[[res$e_col]], na.rm = TRUE) else 0
      sprintf("Median split | n = %d, events = %d", nrow(res$d_km), n_ev)
    })

    output$vb_cox_val <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res) || !res$cox_testable || is.na(res$hr_sd)) {
        "Gated / Untestable"
      } else {
        sprintf("HR = %.2f (p = %.3g)", res$hr_sd, res$cs$P)
      }
    })
    output$vb_cox_sub <- shiny::renderText({
      res <- reactive_res()
      if (is.null(res) || !res$cox_testable || is.na(res$hr_sd)) return("")
      sprintf("95%% CI: [%.2f, %.2f] | C-index = %.2f", res$hr_l_sd, res$hr_u_sd, res$cs$C_index)
    })

    # 5. Violin Plot
    output$plot_violin <- shiny::renderPlot({
      res <- reactive_res()
      if (is.null(res) || !res$w_testable || nrow(res$d_w) == 0) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
              label = "Not enough data for 3-year outcome violin plot\nunder current quality gates (minimum group size requirement not met).",
              colour = "grey40", size = 4.5, fontface = "italic", hjust = 0.5) +
            ggplot2::theme_void() +
            ggplot2::labs(title = if (!is.null(res)) paste0("Outcome Distribution: ", res$score_label) else "Outcome Distribution")
        )
      }

      d_w <- res$d_w
      endpoint <- res$endpoint
      d_w$Outcome <- factor(d_w[[res$group_col]],
        levels = if (endpoint == "OS") c("Dead_3yr", "Alive_3yr") else c("Recurred", "Recurrence-Free")
      )
      x_labels <- if (endpoint == "OS") {
        c("Dead_3yr" = "Dead (≤3 yr)", "Alive_3yr" = "Alive (>3 yr)")
      } else {
        c("Recurred" = "Recurred (≤3 yr)", "Recurrence-Free" = "Recurrence-Free")
      }
      n_tbl <- table(d_w$Outcome)
      n_labels <- paste0(x_labels[names(n_tbl)], "\n(n=", n_tbl, ")")
      names(n_labels) <- names(n_tbl)

      r_text <- if (!is.na(res$ws$r)) sprintf("r = %.2f", res$ws$r) else "r = NA"
      p_text <- if (!is.na(res$ws$p)) sprintf("p = %.3g", res$ws$p) else "p = NA"
      stars <- dtpsig:::.sig_stars(res$ws$p, "ns")

      ggplot2::ggplot(d_w, ggplot2::aes(x = Outcome, y = live_score, fill = Outcome)) +
        ggplot2::geom_violin(alpha = 0.85, trim = FALSE, linewidth = 0.4, colour = "grey30") +
        ggplot2::geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = dtpsig::pub_palette, guide = "none") +
        ggplot2::scale_x_discrete(labels = n_labels) +
        ggplot2::labs(
          title = paste0("3-Year Outcome Distribution"),
          subtitle = sprintf("%s | %s\nWilcoxon %s, %s (%s)", res$cohort_display, res$score_label, r_text, p_text, stars),
          x = "3-Year Outcome",
          y = "Score Value"
        ) +
        dtpsig::pub_theme
    })

    # 6. Kaplan-Meier Plot
    output$plot_km <- shiny::renderPlot({
      res <- reactive_res()
      if (is.null(res) || !res$km_testable || nrow(res$d_km) == 0) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
              label = "Not enough data for Kaplan-Meier survival curve\nunder current quality gates (minimum group size or event count not met).",
              colour = "grey40", size = 4.5, fontface = "italic", hjust = 0.5) +
            ggplot2::theme_void() +
            ggplot2::labs(title = if (!is.null(res)) paste0("Kaplan-Meier: ", res$score_label) else "Kaplan-Meier")
        )
      }

      d_km <- res$d_km
      x_cap <- if (res$endpoint == "OS") res$cfg$os_cutpoint else res$cfg$rfs_cutpoint
      if (is.null(x_cap)) x_cap <- 36

      d_km$time_capped <- pmin(as.numeric(d_km[[res$t_col]]), x_cap)
      d_km$event_val <- as.numeric(d_km[[res$e_col]])
      med_score <- stats::median(d_km$live_score, na.rm = TRUE)
      d_km$Group <- factor(
        ifelse(d_km$live_score > med_score, "High", "Low"),
        levels = c("Low", "High")
      )

      sf <- survival::survfit(survival::Surv(time_capped, event_val) ~ Group, data = d_km)
      strata <- rep(sub(".*=", "", names(sf$strata)), sf$strata)
      cd <- rbind(
        data.frame(time = 0, surv = 1, Group = levels(d_km$Group)),
        data.frame(time = sf$time, surv = sf$surv, Group = strata)
      )
      cd$Group <- factor(cd$Group, levels = c("Low", "High"))

      p_text <- if (!is.na(res$km_p)) sprintf("log-rank p = %.3g (%s)", res$km_p, dtpsig:::.sig_stars(res$km_p, "ns")) else "log-rank p = NA"

      ggplot2::ggplot(cd, ggplot2::aes(x = time, y = surv, colour = Group)) +
        ggplot2::geom_step(linewidth = 0.8) +
        ggplot2::scale_colour_manual(
          values = c(Low = "#4DBBD5", High = "#E64B35"),
          name = "Score Group",
          labels = c(Low = "Low (≤ median)", High = "High (> median)")
        ) +
        ggplot2::scale_x_continuous(
          limits = c(0, x_cap),
          breaks = seq(0, x_cap, by = 12),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          labels = scales::percent_format(accuracy = 1),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          title = "Kaplan-Meier Survival Curves (Median Split)",
          subtitle = sprintf("%s | %s\n%s", res$cohort_display, res$score_label, p_text),
          x = "Months from diagnosis",
          y = "Survival probability"
        ) +
        dtpsig::pub_theme
    })

    # 7. DataTable Output
    output$stats_table <- DT::renderDataTable({
      res <- reactive_res()
      if (is.null(res)) {
        return(DT::datatable(data.frame(Message = "No data available"), rownames = FALSE))
      }

      n_km_ev <- if (nrow(res$d_km) > 0) sum(res$d_km[[res$e_col]], na.rm = TRUE) else 0

      df_stats <- data.frame(
        Test = c(
          "Wilcoxon Rank-Sum (3-yr outcome)",
          "Kaplan-Meier Log-Rank (median split)",
          "Univariate Cox Proportional Hazards (per 1 SD)"
        ),
        Statistic = c(
          if (res$w_testable && !is.na(res$ws$r)) sprintf("r = %.3f", res$ws$r) else "Untestable / Gated",
          if (res$km_testable) "Chi-square (df=1)" else "Untestable / Gated",
          if (res$cox_testable && !is.na(res$hr_sd)) sprintf("HR = %.2f (95%% CI: %.2f - %.2f)", res$hr_sd, res$hr_l_sd, res$hr_u_sd) else "Untestable / Gated"
        ),
        N = c(
          if (res$w_testable) as.character(res$ws$n) else sprintf("%d", nrow(res$d_w)),
          sprintf("%d (events: %d)", nrow(res$d_km), n_km_ev),
          if (res$cox_testable) sprintf("%d (events: %d)", res$cs$N, res$cs$N_events) else sprintf("%d (events: %d)", nrow(res$d_km), n_km_ev)
        ),
        P = c(
          if (res$w_testable && !is.na(res$ws$p)) sprintf("%.3g", res$ws$p) else "NA",
          if (res$km_testable && !is.na(res$km_p)) sprintf("%.3g", res$km_p) else "NA",
          if (res$cox_testable && !is.na(res$cs$P)) sprintf("%.3g", res$cs$P) else "NA"
        ),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        df_stats,
        options = list(
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          ordering = FALSE
        ),
        rownames = FALSE,
        caption = htmltools::tags$caption(
          style = "caption-side: bottom; text-align: left; font-style: italic; color: #6c757d; margin-top: 8px;",
          "Live query - not FDR-corrected; see the static CRC report for BH-adjusted results across the full signature panel."
        )
      )
    })
  })
}
