# =============================================================================
# Script 02: Outcome GWAS Harmonization
# =============================================================================
#
# Downloads FinnGen R12 outcome GWAS summary statistics and harmonizes
# instrument SNPs against each outcome.
#
# Input:  Instruments (data/interim/instruments_*.csv)
# Output: Harmonized datasets per target x outcome (data/processed/)
#
# Pre-registration: Sections 7, 8
# =============================================================================

source("scripts/utils.R")

# =============================================================================
# 1. Load instruments
# =============================================================================

instruments <- fread(file.path(dir_interim, "instruments_all_targets.csv"))
log_msg("Loaded ", nrow(instruments), " instruments across ",
        length(unique(instruments$target)), " targets")

# =============================================================================
# 2. Download / load FinnGen R12 outcome GWAS
# =============================================================================
#
# FinnGen R12 summary statistics are available at:
#   https://www.finngen.fi/en/access_results
#   https://storage.googleapis.com/finngen-public-data-r12/summary_stats/
#
# Format: finngen_R12_{ENDPOINT}.gz
#
# MANUAL STEP: Download the 7 outcome GWAS files and place in data/raw/
# Alternatively, this script will attempt to download via URL.

FINNGEN_BASE_URL <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_"

download_finngen <- function(endpoint_id) {
  filename <- paste0("finngen_R12_", endpoint_id, ".gz")
  filepath <- file.path(dir_raw, filename)

  if (file.exists(filepath)) {
    log_msg("  Already downloaded: ", filename)
    return(filepath)
  }

  url <- paste0(FINNGEN_BASE_URL, endpoint_id, ".gz")
  log_msg("  Downloading: ", url)

  tryCatch({
    download.file(url, filepath, mode = "wb", quiet = TRUE)
    log_gwas_download(
      gwas_id     = paste0("FinnGen_R12_", endpoint_id),
      source      = url,
      filename    = filename,
      description = paste0("FinnGen R12 endpoint: ", endpoint_id)
    )
    filepath
  }, error = function(e) {
    warning("Failed to download ", endpoint_id, ": ", e$message,
            "\nPlease download manually from: ", url)
    NULL
  })
}

# Download all outcomes
outcome_files <- list()
for (i in seq_len(nrow(OUTCOMES_FINNGEN))) {
  eid <- OUTCOMES_FINNGEN$id[i]
  log_msg("Fetching FinnGen R12: ", eid)
  outcome_files[[eid]] <- download_finngen(eid)
}

# =============================================================================
# 3. Format instruments as TwoSampleMR exposure data
# =============================================================================

format_exposure <- function(inst, target_name) {
  dat <- data.frame(
    SNP             = inst$SNP,
    beta.exposure   = inst$BETA,
    se.exposure     = inst$SE,
    pval.exposure   = inst$PVAL,
    effect_allele.exposure = inst$A1,
    other_allele.exposure  = inst$A2,
    eaf.exposure    = inst$EAF,
    exposure        = paste0(target_name, "_LDL"),
    id.exposure     = target_name,
    stringsAsFactors = FALSE
  )
  dat
}

# =============================================================================
# 4. Extract outcome data and harmonize per target x outcome
# =============================================================================

harmonized_list <- list()

for (tgt in unique(instruments$target)) {
  tgt_inst <- instruments[target == tgt]
  exposure_dat <- format_exposure(tgt_inst, tgt)

  for (i in seq_len(nrow(OUTCOMES_FINNGEN))) {
    oid   <- OUTCOMES_FINNGEN$id[i]
    label <- OUTCOMES_FINNGEN$label[i]
    fpath <- outcome_files[[oid]]

    if (is.null(fpath) || !file.exists(fpath)) {
      warning("Skipping ", oid, " for ", tgt, " — file not available")
      next
    }

    log_msg("Harmonizing: ", tgt, " x ", oid)

    # Load outcome GWAS
    outcome_gwas <- fread(fpath)

    # FinnGen R12 expected columns:
    #   rsids, chrom, pos, ref, alt, af_alt, beta, sebeta, pval, ...
    # Extract outcome data for instrument SNPs
    outcome_subset <- outcome_gwas[rsids %in% tgt_inst$SNP]

    if (nrow(outcome_subset) == 0) {
      log_msg("  WARNING: No instrument SNPs found in outcome GWAS. ",
              "Will attempt LD proxy lookup.")
    }

    # Format as TwoSampleMR outcome
    outcome_dat <- data.frame(
      SNP             = outcome_subset$rsids,
      beta.outcome    = outcome_subset$beta,
      se.outcome      = outcome_subset$sebeta,
      pval.outcome    = outcome_subset$pval,
      effect_allele.outcome = outcome_subset$alt,
      other_allele.outcome  = outcome_subset$ref,
      eaf.outcome     = outcome_subset$af_alt,
      outcome         = label,
      id.outcome      = oid,
      stringsAsFactors = FALSE
    )

    # --- LD proxy lookup for missing SNPs (Section 6.5) ---
    missing_snps <- setdiff(tgt_inst$SNP, outcome_subset$rsids)
    if (length(missing_snps) > 0) {
      log_msg("  ", length(missing_snps), " SNPs missing from outcome; ",
              "attempting LD proxy lookup (r2 >= ", LD_PROXY_R2, ")")

      for (ms in missing_snps) {
        proxy <- tryCatch({
          ieugwasr::ld_matrix(ms, pop = "EUR")
          # Simplified — in practice use ld_reflookup or LDlinkR
          NULL
        }, error = function(e) NULL)

        # Proxy lookup is best done via LDlinkR::LDproxy or the
        # TwoSampleMR proxy mechanism. For now, flag and move on.
      }
      log_msg("  Proxy lookup: implement via LDlinkR if needed")
    }

    # --- Harmonize (Section 8) ---
    if (nrow(outcome_dat) > 0) {
      harmonized <- tryCatch({
        TwoSampleMR::harmonise_data(
          exposure_dat = exposure_dat,
          outcome_dat  = outcome_dat,
          action       = 2  # Exclude ambiguous palindromic SNPs (EAF 0.42-0.58)
        )
      }, error = function(e) {
        warning("Harmonization failed for ", tgt, " x ", oid, ": ", e$message)
        NULL
      })

      if (!is.null(harmonized) && nrow(harmonized[harmonized$mr_keep, ]) > 0) {
        harmonized$target   <- tgt
        harmonized$category <- OUTCOMES_FINNGEN$category[i]

        pair_id <- paste0(tolower(tgt), "_", tolower(oid))
        harmonized_list[[pair_id]] <- harmonized

        # Save individual harmonized dataset
        out_file <- file.path(dir_processed,
                              paste0("harmonized_", pair_id, ".csv"))
        fwrite(harmonized, out_file)
        log_msg("  Saved: ", out_file,
                " (", sum(harmonized$mr_keep), " SNPs retained)")
      } else {
        warning("No SNPs retained after harmonization for ", tgt, " x ", oid)
      }
    }
  }
}

# =============================================================================
# 5. Summary
# =============================================================================

log_msg("=== HARMONIZATION SUMMARY ===")
for (pair_id in names(harmonized_list)) {
  h <- harmonized_list[[pair_id]]
  log_msg(sprintf("  %s: %d SNPs harmonized (%d kept, %d palindromic excluded)",
                  pair_id,
                  nrow(h),
                  sum(h$mr_keep),
                  sum(!h$mr_keep & h$palindromic)))
}

# Save combined harmonized data
all_harmonized <- rbindlist(harmonized_list, fill = TRUE)
fwrite(all_harmonized, file.path(dir_processed, "harmonized_all.csv"))
log_msg("Saved: ", file.path(dir_processed, "harmonized_all.csv"))

log_msg("Script 02 complete.")
