# =============================================================================
# Script 06: Power Calculations
# =============================================================================
#
# Computes minimum detectable ORs and required sample sizes per
# target x outcome using the Burgess (2014) formula.
#
# Input:  Instruments (data/interim/), outcome case counts
# Output: Power table (outputs/tables/)
#
# Pre-registration: Section 14
# =============================================================================

source("scripts/utils.R")

# =============================================================================
# 1. Load instruments and compute variance explained
# =============================================================================

instruments <- fread(file.path(dir_interim, "instruments_all_targets.csv"))

# GLGC sample size
GLGC_N <- 1600000

# Compute per-SNP R² and total R² per target
# R² ≈ F / (F + n - 2) ≈ (beta/se)² / ((beta/se)² + n) for large n
instruments[, R2_snp := F_STAT / (F_STAT + GLGC_N - 2)]

target_r2 <- instruments[, .(
  n_snps  = .N,
  R2_total = sum(R2_snp),
  mean_F   = mean(F_STAT),
  min_F    = min(F_STAT)
), by = target]

log_msg("=== INSTRUMENT STRENGTH ===")
for (i in seq_len(nrow(target_r2))) {
  log_msg(sprintf("  %s: %d SNPs, R²=%.4f, mean F=%.1f, min F=%.1f",
                  target_r2$target[i], target_r2$n_snps[i],
                  target_r2$R2_total[i], target_r2$mean_F[i],
                  target_r2$min_F[i]))
}

# =============================================================================
# 2. Outcome sample sizes
# =============================================================================

# FinnGen R12 approximate case/control counts (from registration Section 7.2)
outcome_sizes <- data.frame(
  outcome    = c("AB1_OTHER_SEPSIS", "AB1_STREPTO_SEPSIS", "J10_PNEUMOPNEUMO",
                 "J10_INFLUENZA", "AB1_HERPES_SIMPLEX", "AB1_ZOSTER",
                 "ST19_FRACT_FOREA"),
  n_cases    = c(17133, 3239, 1625, 11558, 4677, 7132, 27907),
  n_controls = c(439048, 480316, 415538, 415538, 480316, 480316, 462982),
  stringsAsFactors = FALSE
)
# Total N
outcome_sizes$n_total <- outcome_sizes$n_cases + outcome_sizes$n_controls

# Effective sample size for binary outcomes:
# N_eff = 4 * n_cases * n_controls / (n_cases + n_controls)
outcome_sizes$n_eff <- with(outcome_sizes,
                            4 * n_cases * n_controls / n_total)

# =============================================================================
# 3. Power calculation: Burgess (2014) formula
# =============================================================================
#
# For two-sample MR with binary outcome:
#   Var(beta_MR) ≈ 1 / (N_outcome * R²_instrument * K * (1-K))
#   where K = case proportion
#   Power = Phi(|log(OR)| / sqrt(Var) - z_alpha)
#
# We compute:
#   (a) Minimum detectable OR at 80% power
#   (b) Whether OR = 0.85 is detectable
#   (c) Required sample size for OR = 0.85

power_mr <- function(or, r2, n_outcome, case_prop, alpha = 0.05) {
  log_or <- log(or)
  var_beta <- 1 / (n_outcome * r2 * case_prop * (1 - case_prop))
  se <- sqrt(var_beta)
  z_alpha <- qnorm(1 - alpha / 2)
  power <- pnorm(abs(log_or) / se - z_alpha)
  return(power)
}

min_detectable_or <- function(r2, n_outcome, case_prop, alpha = 0.05,
                              power_target = 0.80) {
  var_beta <- 1 / (n_outcome * r2 * case_prop * (1 - case_prop))
  se <- sqrt(var_beta)
  z_alpha <- qnorm(1 - alpha / 2)
  z_power <- qnorm(power_target)
  min_log_or <- (z_alpha + z_power) * se
  return(exp(min_log_or))
}

required_n <- function(or, r2, case_prop, alpha = 0.05, power_target = 0.80) {
  log_or <- log(or)
  z_alpha <- qnorm(1 - alpha / 2)
  z_power <- qnorm(power_target)
  n <- ((z_alpha + z_power) / log_or)^2 / (r2 * case_prop * (1 - case_prop))
  return(ceiling(n))
}

# =============================================================================
# 4. Compute power table
# =============================================================================

BONFERRONI_ALPHA <- 0.05 / 21
CLINICAL_OR <- 0.85  # 15% risk reduction per 1 SD LDL-C

power_table <- list()

for (i in seq_len(nrow(target_r2))) {
  tgt  <- target_r2$target[i]
  r2   <- target_r2$R2_total[i]

  for (j in seq_len(nrow(outcome_sizes))) {
    oid  <- outcome_sizes$outcome[j]
    nc   <- outcome_sizes$n_cases[j]
    nt   <- outcome_sizes$n_total[j]

    if (is.na(nc)) next  # skip if case count unknown

    case_prop <- nc / nt

    # Minimum detectable OR at Bonferroni-corrected alpha
    mdo_bonf <- min_detectable_or(r2, nt, case_prop, alpha = BONFERRONI_ALPHA)

    # Minimum detectable OR at nominal alpha
    mdo_nom <- min_detectable_or(r2, nt, case_prop, alpha = 0.05)

    # Power for clinical OR = 0.85 at both thresholds
    pow_bonf <- power_mr(CLINICAL_OR, r2, nt, case_prop,
                         alpha = BONFERRONI_ALPHA)
    pow_nom  <- power_mr(CLINICAL_OR, r2, nt, case_prop, alpha = 0.05)

    # Required N for clinical OR at Bonferroni
    req_n <- required_n(CLINICAL_OR, r2, case_prop, alpha = BONFERRONI_ALPHA)

    # Power-gap flag (Section 14.4)
    underpowered <- mdo_bonf > 1.50 | mdo_bonf < 0.67

    power_table[[paste0(tgt, "_", oid)]] <- data.frame(
      target               = tgt,
      outcome              = oid,
      n_snps               = target_r2$n_snps[i],
      R2                   = r2,
      n_cases              = nc,
      n_total              = nt,
      min_OR_bonferroni    = round(mdo_bonf, 3),
      min_OR_nominal       = round(mdo_nom, 3),
      power_OR085_bonf     = round(pow_bonf, 3),
      power_OR085_nominal  = round(pow_nom, 3),
      required_N_bonf      = req_n,
      underpowered         = underpowered,
      stringsAsFactors     = FALSE
    )
  }
}

power_df <- do.call(rbind, power_table)
rownames(power_df) <- NULL

# =============================================================================
# 5. Summary and export
# =============================================================================

log_msg("=== POWER SUMMARY ===")
for (i in seq_len(nrow(power_df))) {
  flag <- if (power_df$underpowered[i]) " ** UNDERPOWERED **" else ""
  log_msg(sprintf("  %s x %s: min OR = %.3f (Bonf) | power(0.85) = %.1f%%%s",
                  power_df$target[i], power_df$outcome[i],
                  power_df$min_OR_bonferroni[i],
                  power_df$power_OR085_bonf[i] * 100,
                  flag))
}

fwrite(power_df, file.path(dir_tables, "power_calculations.csv"))
log_msg("Saved: ", file.path(dir_tables, "power_calculations.csv"))

log_msg("Script 06 complete.")
