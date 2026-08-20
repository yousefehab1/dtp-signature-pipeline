# dtpsig

Clean, independent rewrite of the CRC "Drug-Tolerant Persister" (DTP) signature
pipeline from my MSc thesis. Same underlying science as the original project
(ssGSEA scoring, survival/DE statistics, gating thresholds, endpoints) —
rebuilt as an installable R package so the static report and the Shiny
dashboard call the exact same functions, with a data-driven signature panel
and a single Excel results workbook in place of the original's scattered
CSVs and duplicated figures.

## Layout

- `R/` — package source: config, caching, ID conversion, signature loading,
  expression prep, ssGSEA scoring, clinical endpoints, statistics, subtyping,
  plotting, Excel export, results-bundle export, analysis modules, and the
  `R/composites/` figure builders.
- `inst/signatures/panel.csv` — **canonical signature panel**. Add or remove a
  gene signature by editing this file (long format: one row per gene per
  signature) — no code change required. `inst/signatures/composite_defs.csv`
  defines which signature pairs form a "Composite" score.
- `data-raw/` — one-off scripts: migrating the old project's signature file
  into the canonical panel, building the symbol→Ensembl mapping, and the
  MSigDB provenance-verification script.
- `pipeline/build_results_bundle.R` — orchestrator that runs every analysis
  module and writes figures, the Excel workbook, and the dashboard's results
  bundle to `output/`.
- `pipeline/pipeline_config.yml` — run-level configuration (thresholds, paths).
  `legacy_cache_dir`/`legacy_gdcdata_dir` are optional: if set to a directory
  with previously-downloaded GEO/GDC data, the pipeline warms its own cache
  from it instead of re-fetching; leave them unset (`null`) on a fresh clone
  and the pipeline downloads what it needs on first run.
- `dashboard/` — Shiny app for interactively exploring results, including a
  live signature multi-select and a custom Composite picker.

## Running

```r
devtools::load_all(".")
source("pipeline/build_results_bundle.R")   # writes output/{figures,DTP_Results.xlsx,results_bundle.rds}

shiny::runApp("dashboard")                  # explore the latest run's results_bundle.rds
```

## Scope note on "dynamic" signatures

ssGSEA cannot be re-run live in a browser session against multi-gigabyte
expression matrices. What the dashboard *does* let you do fully live: filter
which already-scored signatures are shown, and build a custom Composite from
any two of them, with statistics recomputed on the spot. Adding a genuinely
new signature (one the pipeline has never scored) means adding rows to
`inst/signatures/panel.csv` and re-running `pipeline/build_results_bundle.R`
once — no code change, but not instantaneous inside the dashboard.
