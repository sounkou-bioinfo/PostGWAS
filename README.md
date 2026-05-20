# PostGWAS

African-ancestry-first post-GWAS utilities with native PLINK2 `.pgen` support.

Current focus:

- derive population and African-ancestry keep files from `.psam` metadata;
- compute population-weighted LD directly from `.pgen/.pvar/.psam` via `pgenlibr`;
- perform mixed-ancestry summary-statistic imputation without GAUSS-specific intermediate genotype files.
