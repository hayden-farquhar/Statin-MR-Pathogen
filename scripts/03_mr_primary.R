# =============================================================================
# Script 03: Primary MR Analysis — IVW + Core Methods
# =============================================================================
#
# Runs primary MR estimation (IVW, MR-Egger, weighted median, weighted mode)
# for all target x outcome pairs, with heterogeneity statistics and
# multiple testing correction.
#
# Input:  Harmonized datasets (data/processed/harmonized_*.csv)
# Output: Primary results table, forest plots (outputs/)
#
# Pre-registration: Sections 9, 12, 13
# =============================================================================

source("scripts/utils.R")

library(TwoSampleMR)
library(ggplot2)

# =============================================================================
# 1. Load all harmonized data
# =============================================================================

all_harmonized <- fread(file.path(dir_processed, "harmonized_all.csv"))
log_msg("Loaded harmonized data: ", nrow(all_harmonized), " rows across ",
        length(unique(paste(all_harmonized$id.exposure,
                            all_harmonized$id.outcome))), " pairs")

# =============================================================================
# 2. Run MR per target x outcome (Section 9.1)
# =============================================================================
#
# Method selection follows minimum instrument count rules (Section 6.6):
#   n = 0 -> skip
#   n = 1 -> Wald ratio only
#   n = 2 -> IVW only
#   n >= 3 -> full panel

mr_results_list <- list()
het_results_list <- list()

pairs <- unique(all_harmonized[, .(id.exposure, id.outcome, exposure, outcome,
                                    target, category)])

for (j in seq_len(nrow(pairs))) {
  tgt <- pairs$target[j]
  oid <- pairs$id.outcome[j]
  label <- pairs$outcome[j]
  cat_label <- pairs$category[j]

  pair_dat <- all_harmonized[id.exposure == pairs$id.exposure[j] &
                             id.outcome == oid &
                             mr_keep == TRUE]

  n_snps <- nrow(pair_dat)
  log_msg(sprintf("MR: %s x %s (%s) — %d SNPs", tgt, oid, cat_label, n_snps))

  if (n_snps == 0) {
    log_msg("  Skipping: no valid instruments")
    next
  }

  # Select methods based on instrument count
  if (n_snps == 1) {
    methods <- c("mr_wald_ratio")
  } else if (n_snps == 2) {
    methods <- c("mr_ivw")
  } else {
    methods <- c("mr_ivw", "mr_egger_regression",
                 "mr_weighted_median", "mr_weighted_mode")
  }

  # Run MR
  mr_res <- TwoSampleMR::mr(pair_dat, method_list = methods)
  mr_res$target   <- tgt
  mr_res$category <- cat_label
  mr_res$n_snps   <- n_snps
  mr_results_list[[paste0(tgt, "_", oid)]] <- mr_res

  # Heterogeneity (Section 9.2) — only meaningful with n >= 2
  if (n_snps >= 2) {
    het_res <- TwoSampleMR::mr_heterogeneity(pair_dat)
    het_res$target   <- tgt
    het_res$category <- cat_label
    het_results_list[[paste0(tgt, "_", oid)]] <- het_res
  }
}

# =============================================================================
# 3. Combine results
# =============================================================================

mr_results <- rbindlist(mr_results_list, fill = TRUE)
het_results <- if (length(het_results_list) > 0) {
  rbindlist(het_results_list, fill = TRUE)
} else {
  data.table()
}

# =============================================================================
# 4. Multiple testing correction (Section 12)
# =============================================================================
#
# Apply to the 21 primary IVW P-values (3 targets x 7 outcomes)

ivw_results <- mr_results[method %in% c("Inverse variance weighted",
                                          "Wald ratio")]

# FDR (Benjamini-Hochberg)
ivw_results[, fdr_q := p.adjust(pval, method = "BH")]

# Bonferroni
n_tests <- nrow(ivw_results)
bonferroni_alpha <- 0.05 / 21  # pre-registered: 0.05/21
ivw_results[, bonferroni_p := pmin(pval * 21, 1)]
ivw_results[, sig_fdr := fdr_q < 0.05]
ivw_results[, sig_bonferroni := pval < bonferroni_alpha]

# Merge corrections back
mr_results <- merge(mr_results,
                    ivw_results[, .(id.exposure, id.outcome, method,
                                    fdr_q, bonferroni_p, sig_fdr,
                                    sig_bonferroni)],
                    by = c("id.exposure", "id.outcome", "method"),
                    all.x = TRUE)

# =============================================================================
# 5. Pathogen-class heterogeneity testing (Section 13)
# =============================================================================

pathogen_het <- list()

for (tgt in unique(ivw_results$target)) {
  tgt_ivw <- ivw_results[target == tgt & category != "negative_control"]

  if (nrow(tgt_ivw) < 2) next

  # 5a. Cochran's Q across all 6 infection outcomes (Section 13.2.1)
  beta_vals <- tgt_ivw$b
  se_vals   <- tgt_ivw$se
  weights   <- 1 / se_vals^2
  beta_pooled <- sum(weights * beta_vals) / sum(weights)
  Q_stat <- sum(weights * (beta_vals - beta_pooled)^2)
  Q_df   <- length(beta_vals) - 1
  Q_pval <- pchisq(Q_stat, Q_df, lower.tail = FALSE)

  # 5b. Bacterial vs viral interaction contrast (Section 13.2.2)
  bact <- tgt_ivw[category == "bacterial"]
  vir  <- tgt_ivw[category == "viral"]

  if (nrow(bact) > 0 && nrow(vir) > 0) {
    # Inverse-variance pooled estimates per group
    w_bact <- 1 / bact$se^2
    pooled_bact <- sum(w_bact * bact$b) / sum(w_bact)
    se_bact <- 1 / sqrt(sum(w_bact))

    w_vir <- 1 / vir$se^2
    pooled_vir <- sum(w_vir * vir$b) / sum(w_vir)
    se_vir <- 1 / sqrt(sum(w_vir))

    # Difference (bacterial - viral)
    diff_beta <- pooled_bact - pooled_vir
    diff_se   <- sqrt(se_bact^2 + se_vir^2)  # assumes independence
    diff_z    <- diff_beta / diff_se
    diff_pval <- 2 * pnorm(abs(diff_z), lower.tail = FALSE)
  } else {
    diff_beta <- diff_se <- diff_z <- diff_pval <- NA
  }

  pathogen_het[[tgt]] <- data.frame(
    target = tgt,
    Q_stat = Q_stat, Q_df = Q_df, Q_pval = Q_pval,
    contrast_beta = diff_beta, contrast_se = diff_se,
    contrast_z = diff_z, contrast_pval = diff_pval,
    stringsAsFactors = FALSE
  )
}

pathogen_het_df <- do.call(rbind, pathogen_het)

# =============================================================================
# 6. Forest plots
# =============================================================================

# Per-target forest plot across outcomes
for (tgt in unique(mr_results$target)) {
  tgt_ivw <- mr_results[target == tgt &
                         method %in% c("Inverse variance weighted",
                                        "Wald ratio")]

  if (nrow(tgt_ivw) == 0) next

  tgt_ivw[, OR := exp(b)]
  tgt_ivw[, OR_lo := exp(b - 1.96 * se)]
  tgt_ivw[, OR_hi := exp(b + 1.96 * se)]

  p <- ggplot(tgt_ivw, aes(x = OR, y = reorder(outcome, OR),
                            colour = category)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("bacterial" = "#E74C3C",
                                    "viral" = "#3498DB",
                                    "negative_control" = "#95A5A6")) +
    labs(title = paste0(tgt, "-proxied LDL-lowering: MR estimates"),
         subtitle = "OR per 1 SD decrease in LDL-C (IVW/Wald ratio)",
         x = "Odds Ratio (95% CI)", y = NULL, colour = "Pathogen class") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(dir_figures, paste0("forest_", tolower(tgt), ".pdf")),
         p, width = 8, height = 5)
  log_msg("Saved forest plot: ", tgt)
}

# =============================================================================
# 7. Export results
# =============================================================================

fwrite(mr_results, file.path(dir_tables, "mr_primary_results.csv"))
fwrite(het_results, file.path(dir_tables, "mr_heterogeneity.csv"))
fwrite(ivw_results, file.path(dir_tables, "mr_ivw_corrected.csv"))

if (!is.null(pathogen_het_df) && nrow(pathogen_het_df) > 0) {
  fwrite(pathogen_het_df, file.path(dir_tables, "pathogen_heterogeneity.csv"))
}

log_msg("Script 03 complete.")
