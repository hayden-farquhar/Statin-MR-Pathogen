# =============================================================================
# Drug-Target MR for Pathogen-Stratified Statin-in-Sepsis — Shared Utilities
# =============================================================================
#
# Shared configuration, paths, and helper functions for all analysis scripts.
# This file is sourced by every numbered script (01–07).
#
# Pre-registration: https://doi.org/10.17605/OSF.IO/GAV8F
# =============================================================================

library(data.table)
library(TwoSampleMR)
library(ieugwasr)

# Reproducibility seed (per pre-registration Section 18)
SEED <- 60
set.seed(SEED)

# --- Project paths -----------------------------------------------------------
# All paths are relative to the repository root.
# Run scripts from the repository root: Rscript scripts/01_extract_instruments.R

proj_root <- getwd()
dir_raw <- file.path(proj_root, "data", "raw")
dir_interim <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
dir_figures <- file.path(proj_root, "outputs", "figures")
dir_tables <- file.path(proj_root, "outputs", "tables")
dir_supplementary <- file.path(proj_root, "outputs", "supplementary")

for (d in c(dir_raw, dir_interim, dir_processed,
            dir_figures, dir_tables, dir_supplementary)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Drug-target cis-regions (GRCh38, per registration Section 6.1) ----------

CIS_REGIONS <- data.frame(
  target    = c("PCSK9", "HMGCR", "NPC1L1"),
  chr       = c(1L, 5L, 7L),
  gene_start = c(55505221, 74632154, 44552134),
  gene_end   = c(55530525, 74657929, 44580914),
  cis_start  = c(55005221, 74132154, 44052134),
  cis_end    = c(56030525, 75157929, 45080914),
  stringsAsFactors = FALSE
)

# --- Outcome definitions (FinnGen R12, per registration Section 7.2) ---------

OUTCOMES_FINNGEN <- data.frame(
  id        = c("AB1_OTHER_SEPSIS", "AB1_STREPTO_SEPSIS", "J10_PNEUMOPNEUMO",
                "J10_INFLUENZA", "AB1_HERPES_SIMPLEX", "AB1_ZOSTER",
                "ST19_FRACT_FOREA"),
  label     = c("Other septicaemia (gram-neg predominant)",
                "Streptococcal septicaemia",
                "Pneumococcal pneumonia",
                "Influenza (all)", "Herpes simplex infections",
                "Herpes zoster", "Forearm fracture (negative control)"),
  icd10     = c("A41", "A40", "J13", "J09,J10,J11", "B00", "B02",
                "S52"),
  category  = c("bacterial", "bacterial", "bacterial",
                "viral", "viral", "viral", "negative_control"),
  stringsAsFactors = FALSE
)

# --- Instrument QC thresholds (per registration Section 6.2) -----------------

P_THRESHOLD      <- 5e-8
MAF_THRESHOLD    <- 0.01
CLUMP_R2         <- 0.001
CLUMP_KB         <- 10000
F_STAT_THRESHOLD <- 10
LD_PROXY_R2      <- 0.8

# --- Helper: compute F-statistic ---------------------------------------------

compute_f_stat <- function(beta, se, n) {
  (beta / se)^2
}

# --- Helper: log with timestamp -----------------------------------------------

log_msg <- function(...) {
  message(sprintf("[%s] %s", Sys.time(), paste0(...)))
}

# --- GWAS manifest logger -----------------------------------------------------

log_gwas_download <- function(gwas_id, source, filename, description,
                              manifest_path = file.path(dir_raw,
                                                        "gwas_manifest.csv")) {
  entry <- data.frame(
    gwas_id     = gwas_id,
    source      = source,
    filename    = filename,
    description = description,
    download_date = as.character(Sys.Date()),
    download_time = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  if (file.exists(manifest_path)) {
    manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
    manifest <- rbind(manifest, entry)
  } else {
    manifest <- entry
  }
  write.csv(manifest, manifest_path, row.names = FALSE)
  log_msg("Logged GWAS download: ", gwas_id, " -> ", filename)
}
