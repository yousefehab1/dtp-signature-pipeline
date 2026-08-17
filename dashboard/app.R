# ==============================================================================
# dashboard/app.R - Main Shiny application entry point for dtpsig.
# ==============================================================================

# Source global environment if running standalone
if (!exists("bundle", envir = .GlobalEnv)) {
  if (file.exists("global.R")) {
    source("global.R")
  } else if (file.exists("dashboard/global.R")) {
    source("dashboard/global.R")
  }
}

ui <- bslib::page_navbar(
  title = "dtpsig Dashboard",
  theme = bslib::bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#3B6EA5"
  ),
  bslib::nav_panel(
    title = "Overview",
    mod_overview_ui("overview")
  ),
  bslib::nav_panel(
    title = "Signature Explorer",
    mod_signature_explorer_ui("explorer")
  ),
  bslib::nav_panel(
    title = "CRC Survival",
    mod_crc_ui("crc")
  ),
  bslib::nav_panel(
    title = "Pan-Cancer Survival",
    mod_pancan_ui("pancan")
  ),
  bslib::nav_panel(
    title = "Metastasis DE/GSEA",
    mod_mets_ui("mets")
  ),
  bslib::nav_panel(
    title = "Tables",
    mod_tables_ui("tables")
  ),
  bslib::nav_panel(
    title = "About",
    mod_about_ui("about")
  )
)

server <- function(input, output, session) {
  mod_overview_server("overview", bundle = bundle)
  mod_signature_explorer_server("explorer", bundle = bundle, is_stale = is_bundle_stale)
  mod_crc_server("crc", bundle = bundle)
  mod_pancan_server("pancan", bundle = bundle)
  mod_mets_server("mets", bundle = bundle)
  mod_tables_server("tables", bundle = bundle)
  mod_about_server("about", bundle = bundle)
}

shiny::shinyApp(ui = ui, server = server)
