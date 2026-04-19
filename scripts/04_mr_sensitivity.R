# =============================================================================
# Script 04: Sensitivity Analyses
# =============================================================================
#
# MR-PRESSO, Radial MR, leave-one-out, Steiger directionality,
# alternative LD thresholds.
#
# Input:  Harmonized datasets (data/processed/harmonized_all.csv)
# Output: Sensitivity results (outputs/tables/, outputs/supplementary/)
#
# Pre-registration: Section 10
# =============================================================================

source("scripts/utils.R")

library(TwoSampleMR)
library(MRPRESSO)
library(RadialMR)

# =============================================================================
# 1. Load harmonized data
# =============================================================================

all_harmonized <- fread(file.path(dir_processed, "harmonized_all.csv"))

pairs <- unique(all_harmonized[, .(id.exposure, id.outcome, exposure, outcome,
                                    target, category)])

# =============================================================================
# 2. MR-PRESSO (Section 10.1)
# =============================================================================
#
# Global test with 10,000 simulations. If global P < 0.05, report
# outlier-corrected estimate alongside uncorrected.

presso_results <- list()

for (j in seq_len(nrow(pairs))) {
  tgt <- pairs$target[j]
  oid <- pairs$id.outcome[j]

  pair_dat <- all_harmonized[id.exposure == pairs$id.exposure[j] &
                             id.outcome == oid & mr_keep == TRUE]

  n_snps <- nrow(pair_dat)
  if (n_snps < 3) {
    log_msg("MR-PRESSO skipped for ", tgt, " x ", oid, " (n=", n_snps, ")")
    next
  }

  log_msg("MR-PRESSO: ", tgt, " x ", oid)

  presso_input <- data.frame(
    bx  = pair_dat$beta.exposure,
    by  = pair_dat$beta.outcome,
    bxse = pair_dat$se.exposure,
    byse = pair_dat$se.outcome
  )

  res <- tryCatch({
    MRPRESSO::mr_presso(
      BetaOutcome   = "by",
      BetaExposure  = "bx",
      SdOutcome     = "byse",
      SdExposure    = "bxse",
      OUTLIERtest   = TRUE,
      DISTORTIONtest = TRUE,
      data          = presso_input,
      NbDistribution = 10000,
      SignifThreshold = 0.05
    )
  }, error = function(e) {
    warning("MR-PRESSO failed for ", tgt, " x ", oid, ": ", e$message)
    NULL
  })

  if (!is.null(res)) {
    # Extract global test
    global_p <- res$`MR-PRESSO results`$`Global Test`$Pvalue
    main_res <- res$`Main MR results`

    presso_results[[paste0(tgt, "_", oid)]] <- data.frame(
      target          = tgt,
      outcome         = oid,
      global_p        = global_p,
      n_outliers      = sum(res$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05,
                            na.rm = TRUE),
      beta_raw        = main_res$`Causal Estimate`[1],
      se_raw          = main_res$Sd[1],
      p_raw           = main_res$`P-value`[1],
      beta_corrected  = ifelse(nrow(main_res) > 1, main_res$`Causal Estimate`[2], NA),
      se_corrected    = ifelse(nrow(main_res) > 1, main_res$Sd[2], NA),
      p_corrected     = ifelse(nrow(main_res) > 1, main_res$`P-value`[2], NA),
      stringsAsFactors = FALSE
    )
  }
}

if (length(presso_results) > 0) {
  presso_df <- do.call(rbind, presso_results)
  fwrite(presso_df, file.path(dir_tables, "mr_presso_results.csv"))
  log_msg("Saved MR-PRESSO results")
}

# =============================================================================
# 3. Radial MR (Section 10.1)
# =============================================================================

radial_results <- list()

for (j in seq_len(nrow(pairs))) {
  tgt <- pairs$target[j]
  oid <- pairs$id.outcome[j]

  pair_dat <- all_harmonized[id.exposure == pairs$id.exposure[j] &
                             id.outcome == oid & mr_keep == TRUE]

  if (nrow(pair_dat) < 3) next

  log_msg("Radial MR: ", tgt, " x ", oid)

  radial_input <- RadialMR::format_radial(
    BXG  = pair_dat$beta.exposure,
    BYG  = pair_dat$beta.outcome,
    seBXG = pair_dat$se.exposure,
    seBYG = pair_dat$se.outcome,
    RSID  = pair_dat$SNP
  )

  res_ivw <- tryCatch({
    RadialMR::ivw_radial(radial_input, alpha = 0.05)
  }, error = function(e) {
    warning("Radial IVW failed for ", tgt, " x ", oid, ": ", e$message)
    NULL
  })

  if (!is.null(res_ivw)) {
    n_outliers <- length(res_ivw$outliers)
    radial_results[[paste0(tgt, "_", oid)]] <- data.frame(
      target      = tgt,
      outcome     = oid,
      method      = "Radial IVW",
      beta        = res_ivw$coef[1],
      se          = res_ivw$coef[2],
      n_outliers  = n_outliers,
      Q_stat      = res_ivw$qstatistic,
      stringsAsFactors = FALSE
    )

    # Save radial plot
    pdf(file.path(dir_supplementary,
                  paste0("radial_", tolower(tgt), "_", tolower(oid), ".pdf")),
        width = 7, height = 6)
    RadialMR::plot_radial(res_ivw, radial_scale = TRUE)
    dev.off()
  }
}

if (length(radial_results) > 0) {
  radial_df <- do.call(rbind, radial_results)
  fwrite(radial_df, file.path(dir_tables, "radial_mr_results.csv"))
  log_msg("Saved Radial MR results")
}

# =============================================================================
# 4. Leave-one-out (Section 10.2)
# =============================================================================

loo_results_list <- list()

for (j in seq_len(nrow(pairs))) {
  tgt <- pairs$target[j]
  oid <- pairs$id.outcome[j]

  pair_dat <- all_harmonized[id.exposure == pairs$id.exposure[j] &
                             id.outcome == oid & mr_keep == TRUE]

  if (nrow(pair_dat) < 3) next

  log_msg("Leave-one-out: ", tgt, " x ", oid)

  loo <- TwoSampleMR::mr_leaveoneout(pair_dat)
  loo$target <- tgt

  loo_results_list[[paste0(tgt, "_", oid)]] <- loo

  # LOO plot
  p <- TwoSampleMR::mr_leaveoneout_plot(loo)
  ggsave(file.path(dir_supplementary,
                   paste0("loo_", tolower(tgt), "_", tolower(oid), ".pdf")),
         p[[1]], width = 8, height = max(4, 0.4 * nrow(pair_dat) + 2))
}

if (length(loo_results_list) > 0) {
  loo_df <- rbindlist(loo_results_list, fill = TRUE)
  fwrite(loo_df, file.path(dir_tables, "leave_one_out_results.csv"))
  log_msg("Saved leave-one-out results")
}

# =============================================================================
# 5. Steiger directionality filtering (Section 10.3)
# =============================================================================
#
# For each SNP: compare R² for exposure vs outcome. Exclude SNPs where
# outcome R² > exposure R² (reverse causation). Re-run IVW.

steiger_results <- list()

for (j in seq_len(nrow(pairs))) {
  tgt <- pairs$target[j]
  oid <- pairs$id.outcome[j]

  pair_dat <- all_harmonized[id.exposure == pairs$id.exposure[j] &
                             id.outcome == oid & mr_keep == TRUE]

  if (nrow(pair_dat) < 2) next

  log_msg("Steiger filtering: ", tgt, " x ", oid)

  steiger <- TwoSampleMR::directionality_test(pair_dat)

  # Steiger filtering: keep only SNPs with correct causal direction
  pair_dat$steiger_keep <- TwoSampleMR::steiger_filtering(pair_dat)$steiger_dir

  pair_filtered <- pair_dat[pair_dat$steiger_keep == TRUE, ]
  n_removed <- nrow(pair_dat) - nrow(pair_filtered)

  log_msg("  Steiger: ", n_removed, " SNPs removed (reverse causation)")

  # Re-run IVW on filtered set
  if (nrow(pair_filtered) >= 2) {
    mr_steiger <- TwoSampleMR::mr(pair_filtered,
                                   method_list = c("mr_ivw"))
    mr_steiger$target <- tgt
    mr_steiger$n_removed_steiger <- n_removed

    steiger_results[[paste0(tgt, "_", oid)]] <- mr_steiger
  }
}

if (length(steiger_results) > 0) {
  steiger_df <- rbindlist(steiger_results, fill = TRUE)
  fwrite(steiger_df, file.path(dir_tables, "steiger_filtered_results.csv"))
  log_msg("Saved Steiger-filtered results")
}

# =============================================================================
# 6. Alternative LD clumping thresholds (Section 10.4)
# =============================================================================
#
# This requires re-clumping the original instruments at r2 < 0.01 and < 0.05.
# Instruments are re-extracted from the pre-clump set saved in Script 01.
# Implementation deferred to after primary results are confirmed.

log_msg("NOTE: Alternative LD thresholds (r2 < 0.01, < 0.05) require ",
        "re-running Script 01 with modified CLUMP_R2. ",
        "Run separately after primary analysis confirms data integrity.")

log_msg("Script 04 complete.")
