# Raw Data Acquisition

This directory should contain the GWAS summary statistics listed below. These files are publicly available but too large to include in the repository.

## Exposure GWAS

**GLGC 2021 LDL-C (European ancestry)**

- File: `LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results` (or `.gz`)
- Source: Global Lipids Genetics Consortium (Graham et al., *Nature*, 2021)
- Download: http://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/
- Select: `LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz`
- Size: ~2.1 GB compressed, ~5 GB uncompressed
- If downloaded as `.gz`, decompress before running Script 01: `gzip -dk LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz`

## Outcome GWASs

**FinnGen Release 12 summary statistics**

- Source: https://www.finngen.fi/en/access_results
- Base URL: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/`
- Download the following 7 files:

| Filename | Endpoint | ICD-10 | Cases |
|----------|----------|--------|-------|
| `finngen_R12_AB1_OTHER_SEPSIS.gz` | Other septicaemia | A41 | 17,133 |
| `finngen_R12_AB1_STREPTO_SEPSIS.gz` | Streptococcal septicaemia | A40 | 3,239 |
| `finngen_R12_J10_PNEUMOPNEUMO.gz` | Pneumococcal pneumonia | J13 | 1,625 |
| `finngen_R12_J10_INFLUENZA.gz` | Influenza | J09-J11 | 11,558 |
| `finngen_R12_AB1_HERPES_SIMPLEX.gz` | Herpes simplex | B00 | 4,677 |
| `finngen_R12_AB1_ZOSTER.gz` | Herpes zoster | B02 | 7,132 |
| `finngen_R12_ST19_FRACT_FOREA.gz` | Forearm fracture | S52 | 27,907 |

Each file is approximately 770 MB.

## LD Reference Panel

**1000 Genomes Phase 3 EUR**

- Required for LD clumping in Script 01
- Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz (~1.5 GB)
- Extract EUR files: `tar -xzf 1kg.v3.tgz EUR.bed EUR.bim EUR.fam`
- Place `EUR.bed`, `EUR.bim`, `EUR.fam` in `data/reference/`

## FinnGen Manifest (optional)

- File: `finngen_R12_manifest.tsv`
- Download: `curl -o finngen_R12_manifest.tsv https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`
- Used by Script 00 for endpoint resolution

## plink

- plink v1.90 must be installed and on PATH
- Download: https://www.cog-genomics.org/plink/
