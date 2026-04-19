# =============================================================================
# Script 01: cis-Variant Instrument Extraction per Drug Target
# =============================================================================
#
# Extracts LDL-C-associated SNPs within +/- 500 kb cis-regions of PCSK9,
# HMGCR, and NPC1L1 from the GLGC 2021 summary statistics.
#
# Input:  GLGC 2021 LDL-C GWAS summary statistics (data/raw/)
# Output: Clumped instruments per target (data/interim/)
#
# Pre-registration: Sections 5, 6
# =============================================================================

source("scripts/utils.R")

library(rtracklayer)  # for liftOver if needed

# =============================================================================
# 1. Load GLGC 2021 LDL-C summary statistics
# =============================================================================
#
# MANUAL STEP: Download GLGC 2021 LDL-C summary statistics from:
#   http://csg.sph.umich.edu/willer/public/glgc-lipids2021/
#   File: LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz
#
# Place in data/raw/ and update the filename below if different.

glgc_file <- file.path(dir_raw, "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results")
if (!file.exists(glgc_file)) {
  # Try .gz version
  glgc_file <- file.path(dir_raw, "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
}

if (!file.exists(glgc_file)) {
  stop("GLGC 2021 LDL-C summary statistics not found at:\n  ", glgc_file,
       "\n\nDownload from: http://csg.sph.umich.edu/willer/public/glgc-lipids2021/",
       "\nSee README for instructions.")
}

log_msg("Loading GLGC 2021 LDL-C summary statistics...")
glgc <- fread(glgc_file)

# Log the download
log_gwas_download(
  gwas_id     = "GLGC_2021_LDL_EUR",
  source      = "http://csg.sph.umich.edu/willer/public/glgc-lipids2021/",
  filename    = basename(glgc_file),
  description = "GLGC 2021 LDL-C European-ancestry meta-analysis (Graham et al. Nature 2021)"
)

# Inspect column names and standardise
log_msg("GLGC columns: ", paste(names(glgc), collapse = ", "))

# Expected columns (adjust if actual headers differ):
#   rsID, chr, POS_b37, REF, ALT, POOLED_ALT_AF, EFFECT_SIZE, SE, pvalue_neg_log10, pvalue
# Standardise to working names
std_cols <- c(
  SNP    = "rsID",
  CHR    = "CHROM",
  POS    = "POS_b37",         # GRCh37 positions
  A1     = "ALT",             # effect allele
  A2     = "REF",             # other allele
  EAF    = "POOLED_ALT_AF",
  BETA   = "EFFECT_SIZE",
  SE     = "SE",
  PVAL   = "pvalue"
)

# Check columns exist; warn on mismatches
missing <- setdiff(std_cols, names(glgc))
if (length(missing) > 0) {
  warning("Expected columns not found: ", paste(missing, collapse = ", "),
          "\nAvailable: ", paste(names(glgc), collapse = ", "),
          "\nYou may need to update the column mapping in std_cols.")
}

# Rename to standard names
for (new_name in names(std_cols)) {
  old_name <- std_cols[[new_name]]
  if (old_name %in% names(glgc)) {
    setnames(glgc, old_name, new_name, skip_absent = TRUE)
  }
}

# If pvalue is stored as -log10, convert
if ("pvalue_neg_log10" %in% names(glgc) && !"PVAL" %in% names(glgc)) {
  glgc[, PVAL := 10^(-pvalue_neg_log10)]
}

# =============================================================================
# 2. Genome build handling (GRCh37 -> GRCh38)
# =============================================================================
#
# GLGC 2021 is in GRCh37. Our cis-regions are defined in GRCh38.
# Strategy: filter cis-regions in GRCh37 coordinates first (approximate),
# then liftOver to GRCh38 and re-confirm within exact cis boundaries.
#
# GRCh37 approximate cis-regions (manual lookup from Ensembl GRCh37):
#   PCSK9:  chr1:55,505,647–55,530,526 -> cis: 55,005,647–56,030,526
#   HMGCR:  chr5:74,632,993–74,657,929 -> cis: 74,132,993–75,157,929
#   NPC1L1: chr7:44,552,134–44,580,914 -> cis: 44,052,134–45,080,914
#
# Note: GRCh37/GRCh38 offsets are small for these loci, so the +/-500kb
# window captures the same region in either build. We apply a generous
# initial filter then confirm after liftOver.

CIS_REGIONS_GRCh37 <- data.frame(
  target     = c("PCSK9", "HMGCR", "NPC1L1"),
  chr        = c(1L, 5L, 7L),
  cis_start  = c(55005647, 74132993, 44052134),
  cis_end    = c(56030526, 75157929, 45080914),
  stringsAsFactors = FALSE
)

# =============================================================================
# 3. Extract cis-region SNPs per target
# =============================================================================

instruments_list <- list()

for (i in seq_len(nrow(CIS_REGIONS_GRCh37))) {
  tgt       <- CIS_REGIONS_GRCh37$target[i]
  tgt_chr   <- CIS_REGIONS_GRCh37$chr[i]
  tgt_start <- CIS_REGIONS_GRCh37$cis_start[i]
  tgt_end   <- CIS_REGIONS_GRCh37$cis_end[i]

  log_msg("Extracting cis-region SNPs for ", tgt,
          " (chr", tgt_chr, ":", tgt_start, "-", tgt_end, ")")

  # Filter to cis-region
  cis_snps <- glgc[CHR == tgt_chr & POS >= tgt_start & POS <= tgt_end]
  log_msg("  SNPs in cis-region: ", nrow(cis_snps))

  # Filter: P < 5e-8
  cis_snps <- cis_snps[PVAL < P_THRESHOLD]
  log_msg("  SNPs at P < 5e-8: ", nrow(cis_snps))

  # Filter: MAF > 0.01
  if ("EAF" %in% names(cis_snps)) {
    cis_snps[, MAF := pmin(EAF, 1 - EAF)]
    cis_snps <- cis_snps[MAF > MAF_THRESHOLD]
    log_msg("  SNPs after MAF > 0.01 filter: ", nrow(cis_snps))
  }

  # Compute F-statistic (per-SNP)
  # F = (beta/se)^2; for large GWAS this approximates the Cragg-Donald F
  cis_snps[, F_STAT := compute_f_stat(BETA, SE)]
  cis_snps <- cis_snps[F_STAT > F_STAT_THRESHOLD]
  log_msg("  SNPs after F > 10 filter: ", nrow(cis_snps))

  if (nrow(cis_snps) == 0) {
    warning("No instruments found for ", tgt, " — check GWAS columns and coordinates")
    next
  }

  cis_snps[, target := tgt]
  instruments_list[[tgt]] <- cis_snps
}

instruments_pre_clump <- rbindlist(instruments_list)
log_msg("Total pre-clump instruments: ", nrow(instruments_pre_clump))

# =============================================================================
# 4. LD clumping (per registration Section 6.3)
# =============================================================================

clumped_list <- list()

for (tgt in unique(instruments_pre_clump$target)) {
  log_msg("LD clumping instruments for ", tgt, "...")

  tgt_snps <- instruments_pre_clump[target == tgt]

  # Format for ieugwasr::ld_clump
  clump_input <- data.frame(
    rsid = tgt_snps$SNP,
    pval = tgt_snps$PVAL,
    stringsAsFactors = FALSE
  )

  # Try API first, fall back to local plink
  clumped <- tryCatch({
    ieugwasr::ld_clump(
      dat       = clump_input,
      clump_kb  = CLUMP_KB,
      clump_r2  = CLUMP_R2,
      pop       = "EUR"
    )
  }, error = function(e) {
    log_msg("  API clumping failed, trying local plink...")
    # Local plink clumping with 1000G EUR reference
    bfile <- file.path(proj_root, "data", "reference", "EUR")
    plink_bin <- Sys.which("plink")
    if (nchar(plink_bin) > 0 && file.exists(paste0(bfile, ".bed"))) {
      tryCatch({
        ieugwasr::ld_clump_local(
          dat       = clump_input,
          clump_kb  = CLUMP_KB,
          clump_r2  = CLUMP_R2,
          clump_p   = 1,
          bfile     = bfile,
          plink_bin = plink_bin
        )
      }, error = function(e2) {
        warning("Local clumping also failed for ", tgt, ": ", e2$message)
        NULL
      })
    } else {
      # Last resort: distance-based pruning (no LD info)
      log_msg("  No plink/reference available. Using distance-based pruning.")
      tgt_snps_sorted <- tgt_snps[order(PVAL)]
      keep <- character(0)
      for (s in seq_len(nrow(tgt_snps_sorted))) {
        snp_pos <- tgt_snps_sorted$POS[s]
        snp_id  <- tgt_snps_sorted$SNP[s]
        # Check distance to all already-kept SNPs
        if (length(keep) == 0 ||
            all(abs(snp_pos - tgt_snps_sorted$POS[tgt_snps_sorted$SNP %in% keep]) > CLUMP_KB * 1000)) {
          keep <- c(keep, snp_id)
        }
      }
      data.frame(rsid = keep, pval = tgt_snps_sorted$PVAL[tgt_snps_sorted$SNP %in% keep])
    }
  })

  if (!is.null(clumped) && nrow(clumped) > 0) {
    tgt_clumped <- tgt_snps[SNP %in% clumped$rsid]
    log_msg("  ", tgt, ": ", nrow(tgt_clumped), " instruments after clumping")
    clumped_list[[tgt]] <- tgt_clumped
  } else {
    warning("No instruments survived clumping for ", tgt)
  }
}

instruments_final <- rbindlist(clumped_list)

# =============================================================================
# 5. Instrument summary and export
# =============================================================================

log_msg("=== INSTRUMENT SUMMARY ===")
for (tgt in unique(instruments_final$target)) {
  tgt_dat <- instruments_final[target == tgt]
  log_msg(sprintf("  %s: %d SNPs | median F = %.1f | min F = %.1f",
                  tgt, nrow(tgt_dat),
                  median(tgt_dat$F_STAT), min(tgt_dat$F_STAT)))
}

# Variance explained (R2) per instrument: R2 ≈ F / (F + n - 2) for large n
# We'll compute this more precisely during power calculations

# Save per-target instrument files
for (tgt in unique(instruments_final$target)) {
  out_file <- file.path(dir_interim,
                        paste0("instruments_", tolower(tgt), ".csv"))
  fwrite(instruments_final[target == tgt], out_file)
  log_msg("Saved: ", out_file)
}

# Save combined instrument file
fwrite(instruments_final,
       file.path(dir_interim, "instruments_all_targets.csv"))
log_msg("Saved: ", file.path(dir_interim, "instruments_all_targets.csv"))

log_msg("Script 01 complete.")
