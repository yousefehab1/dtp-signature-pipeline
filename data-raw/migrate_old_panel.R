# One-off migration: convert the old thesis project's wide, blank-padded
# signature CSV (Data/final.csv) into this package's canonical long-format
# panel (inst/signatures/panel.csv) plus a composite-definitions table
# (inst/signatures/composite_defs.csv).
#
# Run once, from the package root:
#   Rscript data-raw/migrate_old_panel.R
#
# Not part of the runtime pipeline -- after this has been run, the canonical
# panel is the single source of truth and this script is not needed again
# unless re-deriving the panel from the old project's data.

library(dplyr)
library(tidyr)
library(readr)

old_final_csv <- "/Users/elabd/Master's Work/Thesis/Data/final.csv"
out_panel      <- "inst/signatures/panel.csv"
out_composites <- "inst/signatures/composite_defs.csv"

stopifnot(file.exists(old_final_csv))

# Same BOM-safe header handling as the old project's load_signature_panel()
# (core/signatures.R) -- read.csv()'s fileEncoding="UTF-8-BOM" truncates this
# file at its first non-UTF-8 byte, so the BOM is stripped from the header
# string directly instead.
wide <- read.csv(old_final_csv, stringsAsFactors = FALSE, check.names = FALSE)
names(wide) <- sub("^﻿", "", names(wide))

# Provenance metadata transcribed from the thesis Table 2.4 (Section_2_Methods.md
# section 2.2). Source_Gene_Count is the count as published; it is intentionally
# larger than the number of non-blank rows in the old CSV for several sets,
# because that CSV already reflects post-mapping/dedup loss which this
# migration keeps visible instead of hiding.
metadata <- tribble(
  ~Old_Column, ~Signature_ID, ~Signature_Name,                                  ~Direction,  ~Source_Citation,                                  ~Source_Gene_Count, ~Notes,
  "Up",        "Up",          "DTP up-regulated component",                     "Up",        "Atlasi Laboratory (unpublished)",                 313L,               NA_character_,
  "Down",      "Down",        "DTP down-regulated component",                   "Down",      "Atlasi Laboratory (unpublished)",                 398L,               NA_character_,
  "Fetal",     "Fetal",       "Fetal intestinal programme",                     "Reference", "Mustata et al. 2013, Cell Reports; Yui et al. 2018", 302L,             NA_character_,
  "revSC",     "revSC",       "Revival stem cell programme",                    "Reference", "Ayyaz et al. 2019, Nature",                       396L,               NA_character_,
  "RSC",       "RSC",         "Regenerative stem cell programme",               "Reference", "Vasquez et al. 2022, Cell Stem Cell",              233L,               NA_character_,
  "CSC",       "CBC",         "Crypt-base columnar cell programme",             "Reference", "Vasquez et al. 2022, Cell Stem Cell",              366L,               "Legacy panel column was named 'CSC'; renamed to CBC (crypt-base columnar) per thesis Table 2.4.",
  "MYC",       "MYC",         "MYC target module (Hallmark V1+V2)",             "Reference", "Liberzon et al. 2015, MSigDB Hallmark",           240L,               NA_character_,
  "IBD",       "IBD",         "IBD epithelial inflamed-vs-healthy programme",   "Reference", "Smillie et al. 2019, Cell",                       129L,               NA_character_
)

stopifnot(setequal(metadata$Old_Column, names(wide)))

long <- purrr::map_dfr(metadata$Old_Column, function(col) {
  genes <- unique(as.character(wide[[col]]))
  genes <- genes[!is.na(genes) & trimws(genes) != "" & !tolower(genes) %in% c("na", "none")]
  meta <- metadata[metadata$Old_Column == col, ]
  tibble(
    Signature_ID       = meta$Signature_ID,
    Signature_Name     = meta$Signature_Name,
    Gene_ID            = genes,
    Direction          = meta$Direction,
    Source_Citation    = meta$Source_Citation,
    Source_Gene_Count  = meta$Source_Gene_Count,
    Notes              = meta$Notes
  )
})

# Gene_Symbol is a display/QC convenience filled in once here; it is never
# used for scoring (scoring is always done in the Ensembl namespace, per
# cfg$id_type). Falls back to NA if org.Hs.eg.db is unavailable.
if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                keys = unique(long$Gene_ID),
                                column = "SYMBOL", keytype = "ENSEMBL",
                                multiVals = "first")
  long$Gene_Symbol <- unname(sym[long$Gene_ID])
} else {
  long$Gene_Symbol <- NA_character_
}

long <- long %>%
  select(Signature_ID, Signature_Name, Gene_ID, Gene_Symbol, Direction,
         Source_Citation, Source_Gene_Count, Notes) %>%
  arrange(Signature_ID, Gene_ID)

stopifnot(!anyDuplicated(long[c("Signature_ID", "Gene_ID")]))

dir.create(dirname(out_panel), recursive = TRUE, showWarnings = FALSE)
write_csv(long, out_panel, na = "")

counts <- long %>% count(Signature_ID, name = "N_Genes_Scored") %>% arrange(Signature_ID)
message("Migrated panel gene counts (scored / source):")
print(left_join(counts, distinct(metadata, Signature_ID, Source_Gene_Count), by = "Signature_ID"))

composite_defs <- tribble(
  ~Composite_ID, ~Display_Name,                ~Positive_Signature, ~Negative_Signature, ~Is_Default,
  "Composite",   "DTP Composite (Up - Down)",  "Up",                "Down",              TRUE
)
write_csv(composite_defs, out_composites)

message("Wrote ", out_panel, " and ", out_composites)
