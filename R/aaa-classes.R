# Domain objects -------------------------------------------------------------

.null_or_character <- S7::new_union(NULL, S7::class_character)
.null_or_list <- S7::new_union(NULL, S7::class_list)

.scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    sprintf("@%s must be a non-empty character scalar", name)
  }
}

.scalar_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    sprintf("@%s must be a non-missing numeric scalar", name)
  }
}

#' Genomic locus
#'
#' A small S7 value object identifying a genomic interval. It deliberately keeps
#' biological details such as lead variants and locus construction rules outside
#' the class; those belong to workflow-specific functions.
#'
#' @param chrom Chromosome name.
#' @param start Start position, 1-based.
#' @param end End position, inclusive.
#' @param genome_build Optional genome build, e.g. `"GRCh38"`.
#' @param locus_id Optional stable locus identifier.
#' @export
Locus <- S7::new_class(
  "Locus",
  package = "PostGWAS",
  properties = list(
    chrom = S7::class_character,
    start = S7::class_double,
    end = S7::class_double,
    genome_build = .null_or_character,
    locus_id = .null_or_character
  ),
  validator = function(self) {
    errors <- c(
      .scalar_string(self@chrom, "chrom"),
      .scalar_number(self@start, "start"),
      .scalar_number(self@end, "end")
    )
    if (is.numeric(self@start) && is.numeric(self@end) && self@end < self@start) {
      errors <- c(errors, "@end must be greater than or equal to @start")
    }
    if (!is.null(self@genome_build)) errors <- c(errors, .scalar_string(self@genome_build, "genome_build"))
    if (!is.null(self@locus_id)) errors <- c(errors, .scalar_string(self@locus_id, "locus_id"))
    errors
  }
)

#' Harmonized summary statistics
#'
#' Canonical summary-statistic table plus provenance. The table should use a
#' common allele contract: `effect_allele` is the allele for `beta`/`z`, and
#' `other_allele` is the alternate coding allele.
#'
#' @param data Data frame with at least `variant_id`, `effect_allele`, and
#'   `other_allele`.
#' @param resource_id Optional resource identifier.
#' @param genome_build Optional genome build.
#' @param allele_basis Allele basis for `beta`/`z`, usually `"effect_allele"`.
#' @param provenance List of provenance metadata.
#' @param diagnostics List of QC diagnostics.
#' @export
HarmonizedSumstats <- S7::new_class(
  "HarmonizedSumstats",
  package = "PostGWAS",
  properties = list(
    data = S7::class_any,
    resource_id = .null_or_character,
    genome_build = .null_or_character,
    allele_basis = S7::class_character,
    provenance = S7::class_list,
    diagnostics = S7::class_list
  ),
  validator = function(self) {
    errors <- character()
    if (!is.data.frame(self@data)) {
      errors <- c(errors, "@data must be a data frame")
    } else {
      required <- c("variant_id", "effect_allele", "other_allele")
      missing <- setdiff(required, names(self@data))
      if (length(missing) > 0L) {
        errors <- c(errors, sprintf("@data is missing required column(s): %s", paste(missing, collapse = ", ")))
      }
    }
    errors <- c(errors, .scalar_string(self@allele_basis, "allele_basis"))
    errors
  }
)

#' Dataset-specific LD matrix
#'
#' A signed LD correlation matrix aligned to a variant table for one dataset and
#' one locus. This object carries provenance and diagnostics so LD is never a
#' bare matrix in higher-level coloc workflows.
#'
#' @param dataset_id Optional dataset identifier.
#' @param locus A [Locus].
#' @param variants Variant data frame aligned to rows/columns of `R`; must
#'   contain `variant_id`.
#' @param R Numeric signed LD correlation matrix.
#' @param allele_basis Allele basis for `R`.
#' @param genome_build Optional genome build.
#' @param provider_id Optional LD-provider identifier.
#' @param diagnostics List of LD diagnostics.
#' @param provenance List of provenance metadata.
#' @export
DatasetLD <- S7::new_class(
  "DatasetLD",
  package = "PostGWAS",
  properties = list(
    dataset_id = .null_or_character,
    locus = Locus,
    variants = S7::class_any,
    R = S7::class_any,
    allele_basis = S7::class_character,
    genome_build = .null_or_character,
    provider_id = .null_or_character,
    diagnostics = S7::class_list,
    provenance = S7::class_list
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
      if (is.matrix(self@R) && nrow(self@variants) != nrow(self@R)) {
        errors <- c(errors, "@variants must have one row per row/column of @R")
      }
      if (!"variant_id" %in% names(self@variants)) {
        errors <- c(errors, "@variants must contain a variant_id column")
      }
    }
    errors <- c(errors, .scalar_string(self@allele_basis, "allele_basis"))
    errors
  }
)

#' Ancestry or reference-panel weights
#'
#' A named, non-negative weight vector plus method metadata.
#'
#' @param dataset_id Optional dataset identifier.
#' @param weights Named numeric weight vector.
#' @param method Weight-estimation method.
#' @param reference Optional reference resource identifier.
#' @param diagnostics List of diagnostics.
#' @param provenance List of provenance metadata.
#' @export
AncestryWeights <- S7::new_class(
  "AncestryWeights",
  package = "PostGWAS",
  properties = list(
    dataset_id = .null_or_character,
    weights = S7::class_double,
    method = S7::class_character,
    reference = .null_or_character,
    diagnostics = S7::class_list,
    provenance = S7::class_list
  ),
  validator = function(self) {
    errors <- character()
    if (length(self@weights) == 0L || is.null(names(self@weights)) || any(!nzchar(names(self@weights)))) {
      errors <- c(errors, "@weights must be a named numeric vector")
    }
    if (any(is.na(self@weights)) || any(self@weights < 0)) {
      errors <- c(errors, "@weights must be non-missing and non-negative")
    }
    if (sum(self@weights) <= 0) {
      errors <- c(errors, "@weights must have positive total mass")
    }
    errors <- c(errors, .scalar_string(self@method, "method"))
    errors
  }
)

#' Fine-mapping result
#'
#' A normalized container for outputs from tools such as SuSiEx or MsCAVIAR.
#'
#' @param method Fine-mapping method name.
#' @param data Method-specific normalized data frame.
#' @param variants Variant metadata data frame.
#' @param credible_sets List of credible-set records.
#' @param provenance List of provenance metadata.
#' @param diagnostics List of diagnostics.
#' @export
FineMapResult <- S7::new_class(
  "FineMapResult",
  package = "PostGWAS",
  properties = list(
    method = S7::class_character,
    data = S7::class_any,
    variants = S7::class_any,
    credible_sets = S7::class_list,
    provenance = S7::class_list,
    diagnostics = S7::class_list
  ),
  validator = function(self) {
    errors <- .scalar_string(self@method, "method")
    if (!is.data.frame(self@data)) errors <- c(errors, "@data must be a data frame")
    if (!is.data.frame(self@variants)) errors <- c(errors, "@variants must be a data frame")
    errors
  }
)

#' Colocalization input bundle
#'
#' Prepared inputs for a colocalization engine: harmonized summary statistics,
#' dataset-specific LD objects, and the mapping from traits/features to LD.
#'
#' @param locus A [Locus].
#' @param sumstats Named list of [HarmonizedSumstats] objects.
#' @param ld Named list of [DatasetLD] objects.
#' @param dict_sumstatLD Two-column integer matrix mapping sumstat entries to LD
#'   entries, matching the ColocBoost-style convention.
#' @param diagnostics List of diagnostics.
#' @param provenance List of provenance metadata.
#' @export
ColocInputBundle <- S7::new_class(
  "ColocInputBundle",
  package = "PostGWAS",
  properties = list(
    locus = Locus,
    sumstats = S7::class_list,
    ld = S7::class_list,
    dict_sumstatLD = S7::class_any,
    diagnostics = S7::class_list,
    provenance = S7::class_list
  ),
  validator = function(self) {
    errors <- character()
    if (length(self@sumstats) == 0L) errors <- c(errors, "@sumstats must not be empty")
    if (length(self@ld) == 0L) errors <- c(errors, "@ld must not be empty")
    if (!is.matrix(self@dict_sumstatLD)) errors <- c(errors, "@dict_sumstatLD must be a matrix")
    else if (ncol(self@dict_sumstatLD) != 2L) errors <- c(errors, "@dict_sumstatLD must have two columns")
    errors
  }
)

#' Colocalization result
#'
#' Normalized result container returned by coloc engines.
#'
#' @param engine_id Engine identifier.
#' @param result Normalized result data frame.
#' @param raw_result Raw engine output.
#' @param diagnostics List of diagnostics.
#' @param provenance List of provenance metadata.
#' @export
ColocResult <- S7::new_class(
  "ColocResult",
  package = "PostGWAS",
  properties = list(
    engine_id = S7::class_character,
    result = S7::class_any,
    raw_result = S7::class_any,
    diagnostics = S7::class_list,
    provenance = S7::class_list
  ),
  validator = function(self) {
    errors <- .scalar_string(self@engine_id, "engine_id")
    if (!is.data.frame(self@result)) errors <- c(errors, "@result must be a data frame")
    errors
  }
)
