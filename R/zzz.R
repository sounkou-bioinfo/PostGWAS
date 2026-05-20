.register_postgwas_traits <- function() {
  # Trait implementations are runtime registry side effects in s7contract, so
  # they must be installed when the namespace is loaded rather than relying on
  # package build-time evaluation.
  s7contract::impl_trait(SignedCorrelationLD, DenseLDProvider, replace = TRUE)
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  .register_postgwas_traits()
}
