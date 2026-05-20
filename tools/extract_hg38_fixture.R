#!/usr/bin/env Rscript

# Extract a tiny PLINK2 fixture from an HGDP/1KG-style all_hg38 reference.
#
# Required for actual extraction:
#   - all_hg38.pgen or all_hg38.pgen.zst
#   - all_hg38.pvar, all_hg38.pvar.zst, or all_hg38.pvar.gz
#   - all_hg38.psam
#   - plink2 on PATH
#
# The provided Dropbox URL is used for the .pgen.zst when the file is absent.
# Companion .pvar/.psam files must be supplied next to it or via arguments.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}

work_dir <- get_arg("--work-dir", file.path(tempdir(), "postgwas-hg38"))
out_dir <- get_arg("--out-dir", file.path("inst", "extdata", "hg38_fixture"))
pgen <- get_arg("--pgen", file.path(work_dir, "all_hg38.pgen.zst"))
pvar <- get_arg("--pvar", file.path(work_dir, "all_hg38.pvar.zst"))
psam <- get_arg("--psam", file.path(work_dir, "all_hg38.psam"))
chr <- get_arg("--chr", "1")
from_bp <- as.integer(get_arg("--from-bp", "1000000"))
to_bp <- as.integer(get_arg("--to-bp", "2000000"))
max_samples <- as.integer(get_arg("--max-samples", "200"))
url <- get_arg("--pgen-url", "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1")

if (!nzchar(Sys.which("plink2"))) {
  stop("plink2 is required to extract a valid .pgen fixture; install it or provide it on PATH.", call. = FALSE)
}

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pgen)) {
  message("Downloading .pgen.zst from ", url)
  utils::download.file(url, pgen, mode = "wb", quiet = FALSE)
}
if (!file.exists(pvar)) {
  stop("Missing companion .pvar/.pvar.zst: ", pvar, call. = FALSE)
}
if (!file.exists(psam)) {
  stop("Missing companion .psam: ", psam, call. = FALSE)
}

prefix <- file.path(work_dir, "all_hg38")
if (grepl("[.]zst$", pgen)) {
  pgen_unz <- file.path(work_dir, sub("[.]zst$", "", basename(pgen)))
  if (!file.exists(pgen_unz)) {
    system2("zstd", c("-d", "-f", "--rm", "-o", shQuote(pgen_unz), shQuote(pgen)))
  }
  pgen <- pgen_unz
}

# Build a deterministic sample subset from the first max_samples samples.
psam_lines <- readLines(psam, warn = FALSE)
header <- psam_lines[1]
body <- psam_lines[-1]
keep_file <- file.path(work_dir, "fixture.keep")
if (length(body) > max_samples) body <- body[seq_len(max_samples)]
psam_tab <- utils::read.table(text = c(header, body), header = TRUE, comment.char = "", check.names = FALSE)
names(psam_tab) <- sub("^#", "", names(psam_tab))
if ("FID" %in% names(psam_tab)) {
  keep <- psam_tab[, c("FID", "IID"), drop = FALSE]
} else {
  keep <- psam_tab[, "IID", drop = FALSE]
}
utils::write.table(keep, keep_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

out_prefix <- file.path(out_dir, "hg38_chr_fixture")
cmd <- c(
  "--pgen", shQuote(pgen),
  "--pvar", shQuote(pvar),
  "--psam", shQuote(psam),
  "--chr", chr,
  "--from-bp", as.character(from_bp),
  "--to-bp", as.character(to_bp),
  "--keep", shQuote(keep_file),
  "--make-pgen", "vzs",
  "--out", shQuote(out_prefix)
)
status <- system2("plink2", cmd)
if (!identical(status, 0L)) {
  stop("plink2 fixture extraction failed", call. = FALSE)
}

message("Fixture written to: ", normalizePath(out_dir, mustWork = FALSE))
