# Data Dictionary

## Exposure GWAS (GLGC 2021 LDL-C)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| rsID / SNP | character | dbSNP rsID | — |
| CHROM / CHR | integer | Chromosome number | — |
| POS_b37 / POS | integer | Genomic position (GRCh37) | base pairs |
| REF / A2 | character | Reference (non-effect) allele | — |
| ALT / A1 | character | Alternate (effect) allele | — |
| POOLED_ALT_AF / EAF | numeric | Effect allele frequency (European) | proportion |
| EFFECT_SIZE / BETA | numeric | Effect estimate (per allele, on LDL-C) | SD units |
| SE | numeric | Standard error of effect estimate | SD units |
| pvalue / PVAL | numeric | Association p-value | — |
| N | integer | Sample size for this variant | individuals |
| N_studies | integer | Number of contributing studies | — |

## Outcome GWASs (FinnGen R12)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| #chrom | integer | Chromosome number | — |
| pos | integer | Genomic position (GRCh38) | base pairs |
| ref | character | Reference allele | — |
| alt | character | Alternate allele | — |
| rsids | character | dbSNP rsID | — |
| nearest_genes | character | Nearest gene symbol | — |
| pval | numeric | Association p-value | — |
| mlogp | numeric | -log10(p-value) | — |
| beta | numeric | Log-odds ratio (per allele) | log-OR |
| sebeta | numeric | Standard error of beta | log-OR |
| af_alt | numeric | Alternate allele frequency | proportion |
| af_alt_cases | numeric | Alternate allele frequency in cases | proportion |
| af_alt_controls | numeric | Alternate allele frequency in controls | proportion |

## Intermediate Data

### instruments_all_targets.csv (data/interim/)

One row per instrument SNP after LD clumping.

| Variable | Type | Description |
|----------|------|-------------|
| SNP | character | dbSNP rsID |
| CHR | integer | Chromosome |
| POS | integer | Position (GRCh37) |
| A1 | character | Effect allele |
| A2 | character | Other allele |
| EAF | numeric | Effect allele frequency |
| BETA | numeric | Effect on LDL-C (SD units) |
| SE | numeric | Standard error |
| PVAL | character | P-value |
| MAF | numeric | Minor allele frequency |
| F_STAT | numeric | Instrument F-statistic |
| target | character | Drug target (PCSK9, HMGCR, or NPC1L1) |

### harmonized_all.csv (data/processed/)

One row per instrument SNP x outcome pair after TwoSampleMR harmonisation. Contains all columns from TwoSampleMR::harmonise_data() plus:

| Variable | Type | Description |
|----------|------|-------------|
| target | character | Drug target |
| category | character | Outcome category (bacterial, viral, negative_control) |
