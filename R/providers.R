# Simple concrete providers --------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.unique_schema <- function(x) {
  stats::setNames(vapply(x, function(col) paste(class(col), collapse = "/"), character(1)), names(x))
}

.filter_by_locus <- function(data, locus) {
  if ("chrom" %in% names(data)) {
    data <- data[as.character(data$chrom) == as.character(locus@chrom), , drop = FALSE]
  } else if ("chr" %in% names(data)) {
    data <- data[as.character(data$chr) == as.character(locus@chrom), , drop = FALSE]
  }
  pos_col <- if ("pos" %in% names(data)) "pos" else if ("bp" %in% names(data)) "bp" else NULL
  if (!is.null(pos_col)) {
    data <- data[data[[pos_col]] >= locus@start & data[[pos_col]] <= locus@end, , drop = FALSE]
  }
  data
}

#' Data-frame summary-statistic provider
#'
#' A minimal in-memory provider useful for tests, examples, and adapters that
#' have already collected a locus from a table backend.
#'
#' @param data Summary-statistic data frame.
#' @param metadata Provider metadata list.
#' @export
DataFrameSumstatProvider <- S7::new_class(
  "DataFrameSumstatProvider",
  package = "PostGWAS",
  properties = list(
    data = S7::class_any,
    metadata = S7::class_list
  ),
  validator = function(self) {
    if (!is.data.frame(self@data)) "@data must be a data frame"
  }
)

S7::method(query_sumstats, DataFrameSumstatProvider) <- function(provider, locus, traits = NULL,
                                                                 features = NULL, columns = NULL, ...) {
  out <- .filter_by_locus(provider@data, locus)
  if (!is.null(traits) && "trait_id" %in% names(out)) {
    out <- out[out$trait_id %in% traits, , drop = FALSE]
  }
  if (!is.null(features) && "feature_id" %in% names(out)) {
    out <- out[out$feature_id %in% features, , drop = FALSE]
  }
  if (!is.null(columns)) {
    keep <- intersect(columns, names(out))
    out <- out[, keep, drop = FALSE]
  }

  HarmonizedSumstats(
    data = out,
    resource_id = provider@metadata$resource_id %||% NULL,
    genome_build = provider@metadata$genome_build %||% locus@genome_build,
    allele_basis = provider@metadata$allele_basis %||% "effect_allele",
    provenance = list(provider = "DataFrameSumstatProvider"),
    diagnostics = list(n_rows = nrow(out))
  )
}

S7::method(sumstat_schema, DataFrameSumstatProvider) <- function(provider, ...) {
  .unique_schema(provider@data)
}

S7::method(resource_metadata, DataFrameSumstatProvider) <- function(provider, ...) {
  provider@metadata
}

#' Dense LD provider
#'
#' In-memory signed correlation LD provider. This is mainly for tests, cached LD,
#' and small examples; production backends should implement `compute_ld()` for
#' GAUSS, PlinkingDuck, LDGM, or in-sample genotype sources.
#'
#' @param R Numeric signed LD correlation matrix.
#' @param variants Variant metadata data frame aligned to `R`.
#' @param metadata Provider metadata list.
#' @export
DenseLDProvider <- S7::new_class(
  "DenseLDProvider",
  package = "PostGWAS",
  properties = list(
    R = S7::class_any,
    variants = S7::class_any,
    metadata = S7::class_list
  ),
  validator = function(self) {
    errors <- character()
    if (!is.matrix(self@R) || !is.numeric(self@R)) {
      errors <- c(errors, "@R must be a numeric matrix")
    } else if (nrow(self@R) != ncol(self@R)) {
      errors <- c(errors, "@R must be square")
    }
    if (!is.data.frame(self@variants)) {
      errors <- c(errors, "@variants must be a data frame")
    } else {
      if (is.matrix(self@R) && nrow(self@variants) != nrow(self@R)) errors <- c(errors, "@variants must align to @R")
      if (!"variant_id" %in% names(self@variants)) errors <- c(errors, "@variants must contain variant_id")
    }
    errors
  }
)

S7::method(compute_ld, DenseLDProvider) <- function(provider, locus, variants = NULL,
                                                    allele_basis = "effect_allele",
                                                    dataset_id = NULL, ...) {
  R <- provider@R
  var <- provider@variants
  if (!is.null(variants)) {
    if (!"variant_id" %in% names(variants)) {
      stop("`variants` must contain variant_id", call. = FALSE)
    }
    idx <- match(variants$variant_id, var$variant_id)
    keep <- !is.na(idx)
    idx <- idx[keep]
    var <- variants[keep, , drop = FALSE]
    R <- R[idx, idx, drop = FALSE]
  }
  DatasetLD(
    dataset_id = dataset_id %||% provider@metadata$dataset_id %||% NULL,
    locus = locus,
    variants = var,
    R = R,
    allele_basis = allele_basis,
    genome_build = provider@metadata$genome_build %||% locus@genome_build,
    provider_id = provider@metadata$provider_id %||% "dense",
    diagnostics = list(n_variants = nrow(var)),
    provenance = list(provider = "DenseLDProvider")
  )
}

S7::method(resource_metadata, DenseLDProvider) <- function(provider, ...) {
  provider@metadata
}

#' Known ancestry-weight provider
#'
#' Provider for declared or externally estimated ancestry/reference-panel weights.
#'
#' @param weights Named numeric weight vector.
#' @param metadata Provider metadata list.
#' @export
KnownWeightsProvider <- S7::new_class(
  "KnownWeightsProvider",
  package = "PostGWAS",
  properties = list(
    weights = S7::class_double,
    metadata = S7::class_list
  ),
  validator = function(self) {
    AncestryWeights(
      dataset_id = self@metadata$dataset_id %||% NULL,
      weights = self@weights,
      method = self@metadata$method %||% "known",
      reference = self@metadata$reference %||% NULL,
      diagnostics = list(),
      provenance = list()
    )
    character()
  }
)

S7::method(estimate_weights, KnownWeightsProvider) <- function(provider, sumstats, locus = NULL,
                                                               dataset_id = NULL, ...) {
  w <- provider@weights / sum(provider@weights)
  AncestryWeights(
    dataset_id = dataset_id %||% provider@metadata$dataset_id %||% NULL,
    weights = w,
    method = provider@metadata$method %||% "known",
    reference = provider@metadata$reference %||% NULL,
    diagnostics = list(n_weights = length(w)),
    provenance = list(provider = "KnownWeightsProvider")
  )
}

S7::method(resource_metadata, KnownWeightsProvider) <- function(provider, ...) {
  provider@metadata
}
