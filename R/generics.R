# S7 generics ---------------------------------------------------------------

#' Query summary statistics
#'
#' @param provider Object implementing the `SumstatProvider` interface.
#' @param locus A [Locus].
#' @param traits Optional trait IDs to retain.
#' @param features Optional feature IDs to retain.
#' @param columns Optional columns to project.
#' @param ... Backend-specific arguments.
#' @return A [HarmonizedSumstats] object.
#' @export
query_sumstats <- S7::new_generic(
  "query_sumstats",
  "provider",
  function(provider, locus, traits = NULL, features = NULL, columns = NULL, ...) {
    S7::S7_dispatch()
  }
)

#' Inspect a summary-statistic provider schema
#'
#' @param provider Object implementing the `SumstatProvider` interface.
#' @param ... Backend-specific arguments.
#' @return Provider-specific schema metadata.
#' @export
sumstat_schema <- S7::new_generic(
  "sumstat_schema",
  "provider",
  function(provider, ...) S7::S7_dispatch()
)

#' Inspect provider resource metadata
#'
#' @param provider Any resource provider.
#' @param ... Backend-specific arguments.
#' @return A list of metadata.
#' @export
resource_metadata <- S7::new_generic(
  "resource_metadata",
  "provider",
  function(provider, ...) S7::S7_dispatch()
)

#' Compute dataset-specific LD
#'
#' @param provider Object implementing the `LDProvider` interface.
#' @param locus A [Locus].
#' @param variants Variant table to align to. Must contain `variant_id` when
#'   supplied.
#' @param allele_basis Allele basis required for the returned LD.
#' @param dataset_id Dataset identifier for provenance.
#' @param ... Backend-specific arguments.
#' @return A [DatasetLD] object.
#' @export
compute_ld <- S7::new_generic(
  "compute_ld",
  "provider",
  function(provider, locus, variants = NULL, allele_basis = "effect_allele", dataset_id = NULL, ...) {
    S7::S7_dispatch()
  }
)

#' Estimate ancestry weights
#'
#' @param provider Object implementing the `AncestryWeightsProvider` interface.
#' @param sumstats A [HarmonizedSumstats] object or compatible table.
#' @param locus Optional [Locus].
#' @param dataset_id Optional dataset ID.
#' @param ... Backend-specific arguments.
#' @return An [AncestryWeights] object.
#' @export
estimate_weights <- S7::new_generic(
  "estimate_weights",
  "provider",
  function(provider, sumstats, locus = NULL, dataset_id = NULL, ...) {
    S7::S7_dispatch()
  }
)

#' Run a colocalization engine
#'
#' @param engine Object implementing the `ColocEngine` interface.
#' @param bundle A [ColocInputBundle].
#' @param ... Engine-specific arguments.
#' @return A [ColocResult].
#' @export
run_coloc <- S7::new_generic(
  "run_coloc",
  "engine",
  function(engine, bundle, ...) S7::S7_dispatch()
)

#' Inspect colocalization-engine capabilities
#'
#' @param engine A colocalization engine.
#' @param ... Engine-specific arguments.
#' @return A named list.
#' @export
engine_capabilities <- S7::new_generic(
  "engine_capabilities",
  "engine",
  function(engine, ...) S7::S7_dispatch()
)

#' Normalize a raw colocalization result
#'
#' @param engine A colocalization engine.
#' @param raw_result Raw engine output.
#' @param ... Engine-specific arguments.
#' @return A [ColocResult].
#' @export
normalize_result <- S7::new_generic(
  "normalize_result",
  "engine",
  function(engine, raw_result, ...) S7::S7_dispatch()
)

#' Query fine-mapping results
#'
#' @param provider Object implementing the `FineMapResultProvider` interface.
#' @param locus A [Locus].
#' @param traits Optional trait IDs.
#' @param features Optional feature IDs.
#' @param ... Backend-specific arguments.
#' @return A [FineMapResult].
#' @export
query_finemap_results <- S7::new_generic(
  "query_finemap_results",
  "provider",
  function(provider, locus, traits = NULL, features = NULL, ...) S7::S7_dispatch()
)

#' Inspect a fine-mapping provider schema
#'
#' @param provider Object implementing the `FineMapResultProvider` interface.
#' @param ... Backend-specific arguments.
#' @return Provider-specific schema metadata.
#' @export
finemap_schema <- S7::new_generic(
  "finemap_schema",
  "provider",
  function(provider, ...) S7::S7_dispatch()
)
