# =============================================================================
# Script 00: Resolve FinnGen R12 Fracture Negative Control Endpoint
# =============================================================================
#
# Programmatically queries the FinnGen Risteys API to identify the correct
# fracture endpoint phenocode matching the pre-registered ICD-10 S-code
# definition (any fracture).
#
# The pre-registration (Section 7.2) specifies:
#   "Fracture (negative control) — ICD-10 S-codes (any fracture)"
#   "Resolve exact phenocode at data access; ICD-10 S-code definition is binding"
#
# Output: Prints candidate endpoints and updates utils.R OUTCOMES_FINNGEN
#         if a clear match is found.
# =============================================================================

library(httr)
library(jsonlite)

# =============================================================================
# 1. Query FinnGen Risteys API for fracture-related endpoints
# =============================================================================
#
# Risteys provides a public API for endpoint metadata.
# Base URL: https://risteys.finregistry.fi/api/
# Endpoint search: /api/v1/endpoints?search=fracture

cat("=== FinnGen Fracture Endpoint Resolution ===\n\n")

# Method 1: Search Risteys API
cat("--- Method 1: Risteys API search ---\n")

search_risteys <- function(query) {
  url <- paste0("https://risteys.finregistry.fi/api/v1/endpoints?search=",
                URLencode(query))
  resp <- tryCatch({
    httr::GET(url, httr::timeout(30))
  }, error = function(e) {
    message("API request failed: ", e$message)
    NULL
  })

  if (is.null(resp) || httr::status_code(resp) != 200) {
    message("Risteys API unavailable (status: ",
            ifelse(is.null(resp), "connection failed",
                   httr::status_code(resp)), ")")
    return(NULL)
  }

  content <- httr::content(resp, as = "text", encoding = "UTF-8")
  jsonlite::fromJSON(content)
}

fracture_endpoints <- search_risteys("fracture")

if (!is.null(fracture_endpoints) && length(fracture_endpoints) > 0) {
  # Display all fracture-related endpoints
  cat("\nFinnGen endpoints matching 'fracture':\n")
  cat(sprintf("%-35s %-50s %s\n", "PHENOCODE", "NAME", "CASES"))
  cat(paste(rep("-", 100), collapse = ""), "\n")

  if (is.data.frame(fracture_endpoints)) {
    for (i in seq_len(nrow(fracture_endpoints))) {
      cat(sprintf("%-35s %-50s %s\n",
                  fracture_endpoints$endpoint[i],
                  substr(fracture_endpoints$name[i], 1, 50),
                  ifelse(is.null(fracture_endpoints$n_cases[i]), "?",
                         as.character(fracture_endpoints$n_cases[i]))))
    }
  }
}

# Method 2: Check specific candidate phenocodes directly
cat("\n--- Method 2: Direct endpoint lookup ---\n")

candidates <- c(
  "FRACTURES",           # Most likely: broad fracture composite
  "FX_ANY",              # Alternative broad fracture
  "M13_FRACTURES",       # Musculoskeletal chapter
  "INJURIES_FRACTURES",  # Injury chapter
  "ST19_FRACTURES"       # S/T chapter fractures
)

check_endpoint <- function(phenocode) {
  url <- paste0("https://risteys.finregistry.fi/api/v1/endpoints/",
                URLencode(phenocode))
  resp <- tryCatch({
    httr::GET(url, httr::timeout(15))
  }, error = function(e) NULL)

  if (is.null(resp) || httr::status_code(resp) != 200) {
    return(list(exists = FALSE, phenocode = phenocode))
  }

  content <- httr::content(resp, as = "text", encoding = "UTF-8")
  data <- jsonlite::fromJSON(content)
  return(list(
    exists    = TRUE,
    phenocode = phenocode,
    name      = data$name %||% NA,
    n_cases   = data$stats$n_cases %||% NA,
    icd10     = data$ontology$icd10 %||% NA,
    data      = data
  ))
}

cat("\nChecking candidate phenocodes:\n")
for (cand in candidates) {
  result <- check_endpoint(cand)
  if (result$exists) {
    cat(sprintf("  ✓ %s — EXISTS\n", cand))
    cat(sprintf("    Name:   %s\n", result$name))
    cat(sprintf("    Cases:  %s\n", result$n_cases))
    cat(sprintf("    ICD-10: %s\n",
                ifelse(is.null(result$icd10) || length(result$icd10) == 0,
                       "(check Risteys page)",
                       paste(result$icd10, collapse = ", "))))
  } else {
    cat(sprintf("  ✗ %s — not found\n", cand))
  }
}

# =============================================================================
# 2. Method 3: Scrape the FinnGen summary stats listing
# =============================================================================
#
# FinnGen R12 summary stats are listed at a known GCS bucket.
# We can check which fracture-related files exist.

cat("\n--- Method 3: Check GCS bucket for fracture GWAS files ---\n")

check_gwas_file <- function(endpoint_id) {
  url <- paste0("https://storage.googleapis.com/finngen-public-data-r12/",
                "summary_stats/finngen_R12_", endpoint_id, ".gz")
  resp <- tryCatch({
    httr::HEAD(url, httr::timeout(10))
  }, error = function(e) NULL)

  exists <- !is.null(resp) && httr::status_code(resp) == 200
  size_mb <- NA
  if (exists) {
    size <- httr::headers(resp)$`content-length`
    if (!is.null(size)) size_mb <- round(as.numeric(size) / 1e6, 1)
  }
  list(exists = exists, size_mb = size_mb)
}

cat("\nChecking GWAS file availability:\n")
for (cand in candidates) {
  gwas <- check_gwas_file(cand)
  if (gwas$exists) {
    cat(sprintf("  ✓ finngen_R12_%s.gz — EXISTS (%s MB)\n",
                cand, gwas$size_mb))
  } else {
    cat(sprintf("  ✗ finngen_R12_%s.gz — not found\n", cand))
  }
}

# =============================================================================
# 3. Decision and recommendation
# =============================================================================

cat("\n=== RECOMMENDATION ===\n")
cat("
Your pre-registration binds to 'ICD-10 S-codes (any fracture)'.

Decision rules:
1. If FRACTURES exists and its ICD-10 definition covers S02-S92 broadly,
   use FRACTURES. This is the standard FinnGen broad fracture endpoint.

2. If FRACTURES does not exist in R12, check the Risteys search results
   above for the broadest S-code fracture composite.

3. Whatever phenocode you select, verify on the Risteys endpoint page that
   the ICD-10 codes include S-codes (fracture codes). Document the mapping
   in osf/amendment_log.md if the phenocode differs from 'FRACTURES'.

4. Update OUTCOMES_FINNGEN in scripts/utils.R with the confirmed phenocode.
")

cat("\nScript 00 complete.\n")
