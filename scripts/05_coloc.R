# =============================================================================
# Script 05: Colocalization Analysis
# =============================================================================
#
# Bayesian colocalization (coloc.abf) for each drug-target locus x outcome
# pair, using full regional summary statistics.
#
# Input:  Raw GWAS files (data/raw/), cis-region definitions
# Output: Coloc results table, regional plots (outputs/)
#
# Pre-registration: Section 11
# =============================================================================

source("scripts/utils.R")

library(coloc)

# GRCh37 cis-regions (GLGC uses GRCh37 coordinates)
CIS_REGIONS_GRCh37 <- data.frame(
  target     = c("PCSK9", "HMGCR", "NPC1L1"),
  chr        = c(1L, 5L, 7L),
  cis_start  = c(55005647, 74132993, 44052134),
  cis_end    = c(56030526, 75157929, 45080914),
  stringsAsFactors = FALSE
)

# =============================================================================
# 1. Load exposure GWAS (full summary stats for cis-regions)
# =============================================================================

glgc_file <- file.path(dir_raw, "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results")
if (!file.exists(glgc_file)) {
  glgc_file <- file.path(dir_raw, "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
}
if (!file.exists(glgc_file)) {
  stop("GLGC file not found. Run Script 01 first.")
}

glgc <- fread(glgc_file)

# Standardise column names (same mapping as Script 01)
std_cols <- c(SNP = "rsID", CHR = "CHROM", POS = "POS_b37",
              A1 = "ALT", A2 = "REF", EAF = "POOLED_ALT_AF",
              BETA = "EFFECT_SIZE", SE = "SE", PVAL = "pvalue")
for (nn in names(std_cols)) {
  if (std_cols[[nn]] %in% names(glgc)) setnames(glgc, std_cols[[nn]], nn, skip_absent = TRUE)
}
if ("pvalue_neg_log10" %in% names(glgc) && !"PVAL" %in% names(glgc)) {
  glgc[, PVAL := 10^(-pvalue_neg_log10)]
}

# GLGC sample size (approximate European-ancestry; needed for coloc)
GLGC_N <- 1600000

# =============================================================================
# 2. Colocalization priors (Section 11.2)
# =============================================================================

COLOC_PRIORS <- list(
  primary      = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
  sensitivity1 = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
  sensitivity2 = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4)
)

# =============================================================================
# 3. Run coloc per target x outcome
# =============================================================================

coloc_results <- list()

for (i in seq_len(nrow(CIS_REGIONS_GRCh37))) {
  tgt       <- CIS_REGIONS$target[i]
  tgt_chr   <- CIS_REGIONS_GRCh37$chr[i]  # using GRCh37 for GLGC
  tgt_start <- CIS_REGIONS_GRCh37$cis_start[i]
  tgt_end   <- CIS_REGIONS_GRCh37$cis_end[i]

  # Extract exposure data for cis-region
  exp_region <- glgc[CHR == tgt_chr & POS >= tgt_start & POS <= tgt_end]

  if (nrow(exp_region) == 0) {
    warning("No exposure SNPs in cis-region for ", tgt)
    next
  }

  # Deduplicate exposure region by rsID, remove invalid rows
  exp_clean <- exp_region[!is.na(SNP) & SNP != "" & SE > 0 &
                           !is.na(EAF) & EAF > 0 & EAF < 1][!duplicated(SNP)]

  # Case/control counts from FinnGen R12 manifest
  approx_cases <- c(AB1_OTHER_SEPSIS = 17133, AB1_STREPTO_SEPSIS = 3239,
                    J10_PNEUMOPNEUMO = 1625, J10_INFLUENZA = 11558,
                    AB1_HERPES_SIMPLEX = 4677, AB1_ZOSTER = 7132,
                    ST19_FRACT_FOREA = 27907)
  approx_controls <- c(AB1_OTHER_SEPSIS = 439048, AB1_STREPTO_SEPSIS = 480316,
                       J10_PNEUMOPNEUMO = 415538, J10_INFLUENZA = 415538,
                       AB1_HERPES_SIMPLEX = 480316, AB1_ZOSTER = 480316,
                       ST19_FRACT_FOREA = 462982)

  # Iterate over outcomes
  for (j in seq_len(nrow(OUTCOMES_FINNGEN))) {
    oid <- OUTCOMES_FINNGEN$id[j]
    log_msg("Coloc: ", tgt, " x ", oid)

    # Load outcome GWAS
    out_file <- file.path(dir_raw, paste0("finngen_R12_", oid, ".gz"))
    if (!file.exists(out_file)) {
      warning("Outcome GWAS not found: ", out_file)
      next
    }

    out_gwas <- fread(out_file)
    if ("#chrom" %in% names(out_gwas)) setnames(out_gwas, "#chrom", "chrom")

    # Filter to cis-region (FinnGen R12 uses GRCh38 positions)
    tgt_start_38 <- CIS_REGIONS$cis_start[i]
    tgt_end_38   <- CIS_REGIONS$cis_end[i]

    out_region <- out_gwas[chrom == tgt_chr &
                           pos >= tgt_start_38 & pos <= tgt_end_38 &
                           !is.na(rsids) & rsids != "" &
                           sebeta > 0 & !is.na(af_alt) &
                           af_alt > 0 & af_alt < 1][!duplicated(rsids)]

    if (nrow(out_region) == 0) {
      warning("No outcome SNPs in cis-region for ", tgt, " x ", oid)
      next
    }

    n_cases <- approx_cases[oid]
    n_controls <- approx_controls[oid]

    # Merge exposure and outcome by rsID to guarantee alignment
    merged <- merge(
      exp_clean[, .(SNP, POS, BETA, SE, EAF)],
      out_region[, .(rsids, pos, beta, sebeta, af_alt)],
      by.x = "SNP", by.y = "rsids"
    )

    if (nrow(merged) < 50) {
      warning("Only ", nrow(merged), " shared SNPs for ",
              tgt, " x ", oid, " — coloc may be unreliable")
    }
    if (nrow(merged) == 0) next

    log_msg("  ", nrow(merged), " shared SNPs")

    # Use GRCh38 position for both datasets (coloc merges on snp + position)
    exp_coloc_shared <- list(
      snp      = as.character(merged$SNP),
      position = as.integer(merged$pos),
      beta     = as.numeric(merged$BETA),
      varbeta  = as.numeric(merged$SE)^2,
      MAF      = pmin(as.numeric(merged$EAF), 1 - as.numeric(merged$EAF)),
      N        = GLGC_N,
      type     = "quant",
      sdY      = 1
    )

    out_coloc <- list(
      snp      = as.character(merged$SNP),
      position = as.integer(merged$pos),
      beta     = as.numeric(merged$beta),
      varbeta  = as.numeric(merged$sebeta)^2,
      MAF      = pmin(as.numeric(merged$af_alt), 1 - as.numeric(merged$af_alt)),
      N        = n_cases + n_controls,
      s        = n_cases / (n_cases + n_controls),
      type     = "cc"
    )

    # Run coloc with all three prior sets
    for (prior_name in names(COLOC_PRIORS)) {
      priors <- COLOC_PRIORS[[prior_name]]

      res <- tryCatch({
        coloc::coloc.abf(
          dataset1 = exp_coloc_shared,
          dataset2 = out_coloc,
          p1  = priors$p1,
          p2  = priors$p2,
          p12 = priors$p12
        )
      }, error = function(e) {
        warning("Coloc failed for ", tgt, " x ", oid,
                " (", prior_name, "): ", e$message)
        NULL
      })

      if (!is.null(res)) {
        pp <- res$summary
        coloc_results[[paste0(tgt, "_", oid, "_", prior_name)]] <- data.frame(
          target     = tgt,
          outcome    = oid,
          priors     = prior_name,
          n_snps     = pp["nsnps"],
          PP.H0      = pp["PP.H0.abf"],
          PP.H1      = pp["PP.H1.abf"],
          PP.H2      = pp["PP.H2.abf"],
          PP.H3      = pp["PP.H3.abf"],
          PP.H4      = pp["PP.H4.abf"],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# =============================================================================
# 4. Combine and interpret
# =============================================================================

if (length(coloc_results) > 0) {
  coloc_df <- do.call(rbind, coloc_results)
  rownames(coloc_df) <- NULL

  # Interpretation flags (Section 11.3)
  coloc_df$interpretation <- ifelse(
    coloc_df$PP.H4 > 0.8, "Shared causal variant (supports MR)",
    ifelse(coloc_df$PP.H3 > 0.8, "Distinct variants (LD confounding possible)",
           ifelse(coloc_df$PP.H4 > 0.5, "Suggestive (inconclusive)",
                  "No strong evidence"))
  )

  fwrite(coloc_df, file.path(dir_tables, "coloc_results.csv"))
  log_msg("Saved colocalization results")

  # Summary
  log_msg("=== COLOC SUMMARY (primary priors) ===")
  primary_coloc <- coloc_df[coloc_df$priors == "primary", ]
  for (k in seq_len(nrow(primary_coloc))) {
    log_msg(sprintf("  %s x %s: H4=%.3f H3=%.3f — %s",
                    primary_coloc$target[k], primary_coloc$outcome[k],
                    primary_coloc$PP.H4[k], primary_coloc$PP.H3[k],
                    primary_coloc$interpretation[k]))
  }
} else {
  log_msg("No colocalization results generated")
}

log_msg("Script 05 complete.")
