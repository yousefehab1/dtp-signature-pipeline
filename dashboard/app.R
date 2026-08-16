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
  # Phase 5 MVP Tab: Signature Explorer
  bslib::nav_panel(
    title = "Signature Explorer",
    mod_signature_explorer_ui("explorer")
  )
  # Future tabs (Overview, CRC Survival, Pan-Cancer, Mets, Tables, About) will be added here in Phase 7
)

server <- function(input, output, session) {
  # Server logic for Signature Explorer module
  mod_signature_explorer_server("explorer", bundle = bundle, is_stale = is_bundle_stale)
}

shiny::shinyApp(ui = ui, server = server)
