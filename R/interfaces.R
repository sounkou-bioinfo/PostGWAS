# Behavioral contracts -------------------------------------------------------

#' Summary-statistic provider contract
#'
#' `SumstatProvider` is a small behavioral contract for objects that can return
#' harmonized summary statistics for a locus. Concrete implementations may read
#' data frames, files, DuckDB relations, Rduckhts/Rducks tables, or remote stores.
#'
#' @export
SumstatProvider <- s7contract::new_interface(
  "SumstatProvider",
  package = "PostGWAS",
  generics = list(
    query_sumstats = s7contract::interface_requirement(
      query_sumstats,
      args = list(locus = Locus),
      returns = HarmonizedSumstats
    ),
    sumstat_schema = s7contract::interface_requirement(sumstat_schema),
    resource_metadata = s7contract::interface_requirement(resource_metadata)
  )
)

#' LD provider contract
#'
#' `LDProvider` is the core contract for dataset-specific signed LD retrieval.
#' Implementations should return [DatasetLD] objects carrying allele basis,
#' genome build, diagnostics, and provenance.
#'
#' @export
LDProvider <- s7contract::new_interface(
  "LDProvider",
  package = "PostGWAS",
  generics = list(
    compute_ld = s7contract::interface_requirement(
      compute_ld,
      args = list(locus = Locus),
      returns = DatasetLD
    ),
    resource_metadata = s7contract::interface_requirement(resource_metadata)
  )
)

#' Ancestry-weight provider contract
#'
#' `AncestryWeightsProvider` estimates or supplies reference-panel weights for a
#' dataset. GAUSS afmix/zmix, metadata-derived, and manual providers should all
#' satisfy this same small contract.
#'
#' @export
AncestryWeightsProvider <- s7contract::new_interface(
  "AncestryWeightsProvider",
  package = "PostGWAS",
  generics = list(
    estimate_weights = s7contract::interface_requirement(
      estimate_weights,
      returns = AncestryWeights
    ),
    resource_metadata = s7contract::interface_requirement(resource_metadata)
  )
)

#' Colocalization engine contract
#'
#' `ColocEngine` consumes a prepared [ColocInputBundle] and returns a normalized
#' [ColocResult]. ColocBoost, coloc, HyPrColoc, and MACA-style adapters should
#' use this contract rather than exposing unrelated raw call signatures upstream.
#'
#' @export
ColocEngine <- s7contract::new_interface(
  "ColocEngine",
  package = "PostGWAS",
  generics = list(
    run_coloc = s7contract::interface_requirement(
      run_coloc,
      args = list(bundle = ColocInputBundle),
      returns = ColocResult
    ),
    engine_capabilities = s7contract::interface_requirement(engine_capabilities),
    normalize_result = s7contract::interface_requirement(
      normalize_result,
      returns = ColocResult
    )
  )
)

#' Fine-mapping result provider contract
#'
#' `FineMapResultProvider` is for tools whose coloc layer consumes fine-mapping
#' outputs rather than raw summary statistics, such as SuSiEx/MsCAVIAR outputs
#' used by MACA.
#'
#' @export
FineMapResultProvider <- s7contract::new_interface(
  "FineMapResultProvider",
  package = "PostGWAS",
  generics = list(
    query_finemap_results = s7contract::interface_requirement(
      query_finemap_results,
      args = list(locus = Locus),
      returns = FineMapResult
    ),
    finemap_schema = s7contract::interface_requirement(finemap_schema),
    resource_metadata = s7contract::interface_requirement(resource_metadata)
  )
)

# Capability traits ----------------------------------------------------------

#' Signed correlation LD trait
#'
#' Explicit capability marker for LD objects/providers that promise signed LD on
#' the correlation scale. This semantic promise is stronger than merely having a
#' method named `compute_ld()`.
#'
#' @export
SignedCorrelationLD <- s7contract::new_trait(
  "SignedCorrelationLD",
  package = "PostGWAS",
  assoc_consts = list(
    REPRESENTATION = "correlation",
    SIGNED = TRUE
  )
)

#' Sparse precision LD trait
#'
#' Explicit capability marker for LDGM-style sparse precision representations.
#' These objects should be materialized to correlation-scale LD before engines
#' such as ColocBoost consume them.
#'
#' @export
SparsePrecisionLD <- s7contract::new_trait(
  "SparsePrecisionLD",
  package = "PostGWAS",
  assoc_consts = list(
    REPRESENTATION = "precision",
    SIGNED = TRUE
  )
)
