# PostGWAS

PostGWAS is an ancestry-aware post-GWAS workflow package built around small
runtime contracts. Low-level GAUSS-compatible PLINK2 support lives in
`GAUSS.pgen`; PostGWAS adds higher-level resource, LD, weighting, and coloc
abstractions.

Current focus:

- consume `GAUSS.pgen` providers for native `.pgen/.pvar/.psam` LD and
  ancestry-weight estimation;
- define S7/s7contract interfaces for summary statistics, LD providers,
  ancestry-weight providers, coloc engines, and fine-mapping result providers;
- keep LD as a dataset-specific object with allele basis, genome build,
  diagnostics, and provenance rather than passing bare matrices;
- make DuckDB/Rduckhts/Rducks, PlinkingDuck, LDGM, ColocBoost, coloc, HyPrColoc,
  and MACA future adapters instead of hard dependencies.

Core contracts:

- `SumstatProvider`: query harmonized summary statistics for a locus.
- `LDProvider`: return a signed, dataset-specific `DatasetLD` object.
- `AncestryWeightsProvider`: supply or estimate ancestry/reference weights.
- `ColocEngine`: run a coloc method on a prepared `ColocInputBundle`.
- `FineMapResultProvider`: expose SuSiEx/MsCAVIAR-style fine-mapping outputs.

Minimal in-memory adapters are included for development and tests:

- `DataFrameSumstatProvider`
- `DenseLDProvider`
- `KnownWeightsProvider`
