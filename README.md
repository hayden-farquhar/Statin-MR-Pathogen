# Drug-Target Mendelian Randomization Reveals a Power Gap for Pathogen-Stratified Statin-in-Sepsis Decisions

Code repository for: **Drug-Target Mendelian Randomization Reveals a Power Gap for Pathogen-Stratified Statin-in-Sepsis Decisions: A Pre-Registered Analysis**

Hayden Farquhar MBBS MPHTM, Independent Researcher, Finley, NSW, Australia

Pre-registration: [doi:10.17605/OSF.IO/GAV8F](https://doi.org/10.17605/OSF.IO/GAV8F)

Preprint: to be posted

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19648787.svg)](https://doi.org/10.5281/zenodo.19648787)

## Overview

This repository contains the analysis code for a pre-registered drug-target Mendelian randomization study testing whether genetically proxied LDL-lowering via PCSK9, HMGCR, and NPC1L1 has differential causal effects on pathogen-stratified infection outcomes. All 21 analyses were underpowered; the paper's main contribution is a quantitative roadmap of the GWAS sample sizes required to resolve this question.

## Data Sources

| Source | URL | Access | Description |
|--------|-----|--------|-------------|
| GLGC 2021 LDL-C | http://csg.sph.umich.edu/willer/public/glgc-lipids2021/ | Free download | European-ancestry LDL-C GWAS (~1.6M individuals) |
| FinnGen R12 | https://www.finngen.fi/en/access_results | Free download | 7 outcome endpoints (6 infection + 1 negative control) |
| 1000 Genomes EUR | http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz | Free download | LD reference panel for clumping |

Raw data are publicly available but too large (~10 GB total) to include in this repository. See `data/raw/README.md` for download instructions.

## Requirements

- **R** >= 4.3.0
- **plink** v1.90 (on PATH)
- R packages managed via `renv` (versions locked in `renv.lock`)

```bash
# Install R dependencies
R -e 'install.packages("renv"); renv::restore()'
```

## Reproduction

Run all scripts from the repository root directory. Scripts must be run in order — each depends on outputs from preceding scripts.

```bash
# Step 0: Resolve FinnGen fracture endpoint (documents the endpoint selection)
Rscript scripts/00_resolve_fracture_endpoint.R

# Step 1: Extract cis-variant instruments from GLGC 2021
Rscript scripts/01_extract_instruments.R

# Step 2: Download FinnGen R12 outcomes and harmonise
Rscript scripts/02_harmonize.R

# Step 3: Primary MR analysis (IVW, Egger, weighted median/mode)
Rscript scripts/03_mr_primary.R

# Step 4: Sensitivity analyses (MR-PRESSO, radial MR, leave-one-out, Steiger)
Rscript scripts/04_mr_sensitivity.R

# Step 5: Bayesian colocalization (coloc.abf)
Rscript scripts/05_coloc.R

# Step 6: Power calculations and required sample sizes
Rscript scripts/06_power_calc.R

# Step 7: Generate publication figures
Rscript scripts/07_figures.R
```

**Estimated runtime:** ~90 minutes total on a modern laptop (dominated by I/O for loading ~770 MB GWAS files in Scripts 02 and 05). Script 01 requires internet access for LD clumping via the OpenGWAS API (falls back to local plink if unavailable).

## Script Descriptions

| Script | Description | Inputs | Outputs |
|--------|-------------|--------|---------|
| `utils.R` | Shared configuration: paths, cis-regions, outcome definitions, QC thresholds, helper functions | — | — |
| `00_resolve_fracture_endpoint.R` | Programmatic lookup of FinnGen R12 fracture endpoints to resolve pre-registered negative control | FinnGen API / GCS bucket | Console output (decision documentation) |
| `01_extract_instruments.R` | Extract cis-region SNPs for PCSK9, HMGCR, NPC1L1 from GLGC 2021; LD clump; filter on F>10 | `data/raw/LDL_INV_EUR_*.results` | `data/interim/instruments_*.csv` |
| `02_harmonize.R` | Download FinnGen R12 outcome GWASs; harmonise instruments against each outcome | `data/interim/instruments_all_targets.csv`, FinnGen R12 files | `data/processed/harmonized_*.csv` |
| `03_mr_primary.R` | IVW + sensitivity MR methods; FDR/Bonferroni correction; pathogen-class heterogeneity test; forest plots | `data/processed/harmonized_all.csv` | `outputs/tables/mr_primary_results.csv`, `outputs/tables/mr_ivw_corrected.csv`, `outputs/tables/pathogen_heterogeneity.csv`, `outputs/figures/forest_*.pdf` |
| `04_mr_sensitivity.R` | MR-PRESSO, radial MR, leave-one-out, Steiger filtering | `data/processed/harmonized_all.csv` | `outputs/tables/mr_presso_results.csv`, `outputs/tables/radial_mr_results.csv`, `outputs/tables/leave_one_out_results.csv`, `outputs/supplementary/radial_*.pdf` |
| `05_coloc.R` | Bayesian colocalization (coloc.abf) with 3 prior sets per target x outcome | GLGC 2021, FinnGen R12 files | `outputs/tables/coloc_results.csv` |
| `06_power_calc.R` | Burgess (2014) power calculations; minimum detectable ORs; required sample sizes | `data/interim/instruments_all_targets.csv` | `outputs/tables/power_calculations.csv`, `outputs/tables/power_calculations_enhanced.csv` |
| `07_figures.R` | Publication-quality figures: combined forest plot, power-gap bar chart, power simulation curves | `outputs/tables/mr_ivw_corrected.csv`, `outputs/tables/power_calculations_enhanced.csv` | `outputs/figures/figure1_forest_combined.*`, `outputs/figures/figure2_power_gap.*`, `outputs/figures/figure3_power_simulation.*` |

## Outputs

| File | Paper reference |
|------|----------------|
| `outputs/figures/figure1_forest_combined.pdf` | Figure 1 (forest plot) |
| `outputs/figures/figure2_power_gap.pdf` | Figure 2 (required vs current N) |
| `outputs/figures/figure3_power_simulation.pdf` | Figure 3 (power simulation) |
| `outputs/tables/mr_ivw_corrected.csv` | Table 2 (primary MR results) |
| `outputs/tables/power_calculations_enhanced.csv` | Table 3 (required sample sizes) |
| `outputs/tables/mr_primary_results.csv` | Supplementary Table S4 (all MR methods) |
| `outputs/tables/mr_presso_results.csv` | Supplementary Table S5 |
| `outputs/tables/radial_mr_results.csv` | Supplementary Table S6 |
| `outputs/tables/coloc_results.csv` | Supplementary Table S7 |
| `outputs/tables/pathogen_heterogeneity.csv` | Supplementary Table S8 |

## Pre-registration

The full study protocol was pre-registered and frozen on the Open Science Framework prior to data access: [doi:10.17605/OSF.IO/GAV8F](https://doi.org/10.17605/OSF.IO/GAV8F). Three protocol deviations (all identified before viewing results) are documented in the accompanying manuscript.

## Citation

If you use this code, please cite:

```
Farquhar H. Drug-Target Mendelian Randomization Reveals a Power Gap for
Pathogen-Stratified Statin-in-Sepsis Decisions: A Pre-Registered Analysis.
2026. Pre-registration: doi:10.17605/OSF.IO/GAV8F
```

## License

Code: MIT License. Data dictionary and documentation: CC-BY 4.0.
