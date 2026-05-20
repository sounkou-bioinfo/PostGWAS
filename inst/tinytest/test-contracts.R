locus <- Locus(chrom = "1", start = 100, end = 200, genome_build = "GRCh38")

sumstats <- data.frame(
  chrom = c("1", "1", "1"),
  pos = c(90, 150, 210),
  variant_id = c("1:90:A:G", "1:150:A:G", "1:210:A:G"),
  effect_allele = c("G", "G", "G"),
  other_allele = c("A", "A", "A"),
  beta = c(0.1, 0.2, 0.3),
  se = c(0.05, 0.05, 0.05),
  z = c(2, 4, 6),
  trait_id = c("t1", "t1", "t2"),
  stringsAsFactors = FALSE
)

sp <- DataFrameSumstatProvider(
  data = sumstats,
  metadata = list(resource_id = "mock", genome_build = "GRCh38")
)
ht <- query_sumstats(sp, locus)
expect_true(s7contract::implements(sp, SumstatProvider))
expect_true(S7::S7_inherits(ht, HarmonizedSumstats))
expect_equal(nrow(ht@data), 1L)
expect_equal(ht@data$variant_id, "1:150:A:G")

variants <- data.frame(
  variant_id = c("v1", "v2"),
  effect_allele = c("A", "C"),
  other_allele = c("G", "T"),
  stringsAsFactors = FALSE
)
ld_provider <- DenseLDProvider(
  R = matrix(c(1, 0.25, 0.25, 1), 2, 2),
  variants = variants,
  metadata = list(provider_id = "dense-test", genome_build = "GRCh38")
)
ld <- compute_ld(ld_provider, locus, variants = variants, dataset_id = "trait1")
expect_true(s7contract::implements(ld_provider, LDProvider))
expect_true(s7contract::has_trait(ld_provider, SignedCorrelationLD))
expect_equal(dim(ld@R), c(2L, 2L))
expect_equal(ld@dataset_id, "trait1")

wp <- KnownWeightsProvider(
  weights = c(AFR = 8, EUR = 2),
  metadata = list(method = "known", reference = "declared")
)
weights <- estimate_weights(wp, ht, dataset_id = "trait1")
expect_true(s7contract::implements(wp, AncestryWeightsProvider))
expect_equal(sum(weights@weights), 1)
expect_equal(unname(weights@weights["AFR"]), 0.8)
