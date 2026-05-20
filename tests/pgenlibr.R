if (requireNamespace("pgenlibr", quietly = TRUE)) {
  library(PostGWAS)

  extdir <- system.file("extdata", package = "GAUSS.pgen")
  pgen_file <- file.path(extdir, "chr22_hg38_fixture.pgen")
  pvar_file <- file.path(extdir, "chr22_hg38_fixture.pvar.zst")
  psam_file <- file.path(extdir, "chr22_hg38_fixture.psam")

  if (file.exists(pgen_file) && file.exists(pvar_file) && file.exists(psam_file)) {
    pvar <- pgenlibr::NewPvar(pvar_file)
    pgen <- pgenlibr::NewPgen(pgen_file, pvar = pvar)
    on.exit(try(pgenlibr::ClosePgen(pgen), silent = TRUE), add = TRUE)
    on.exit(try(pgenlibr::ClosePvar(pvar), silent = TRUE), add = TRUE)

    nvar <- pgenlibr::GetVariantCt(pvar)
    idx <- seq_len(min(120L, nvar))

    # Orientation invariant: pgenlibr::ReadIntList() is ALT-count dosage.
    # ReadHardcalls(..., allele_num=2) is ALT count; allele_num=1 is REF count.
    orient_idx <- NA_integer_
    for (i in idx) {
      g <- pgenlibr::ReadIntList(pgen, i)[, 1]
      if (length(unique(g[!is.na(g)])) > 1L) {
        orient_idx <- i
        break
      }
    }
    stopifnot(!is.na(orient_idx))
    alt_from_list <- pgenlibr::ReadIntList(pgen, orient_idx)[, 1]
    alt_buf <- pgenlibr::IntBuf(pgen)
    ref_buf <- pgenlibr::IntBuf(pgen)
    pgenlibr::ReadHardcalls(pgen, alt_buf, orient_idx, allele_num = 2L)
    pgenlibr::ReadHardcalls(pgen, ref_buf, orient_idx, allele_num = 1L)
    ok <- !is.na(alt_from_list) & !is.na(alt_buf) & !is.na(ref_buf)
    stopifnot(identical(as.integer(alt_from_list[ok]), as.integer(alt_buf[ok])))
    stopifnot(all(as.integer(alt_buf[ok]) + as.integer(ref_buf[ok]) == 2L))

    sumstats <- data.frame(
      rsid = vapply(idx, function(i) pgenlibr::GetVariantId(pvar, i), character(1)),
      chr = vapply(idx, function(i) pgenlibr::GetVariantChrom(pvar, i), character(1)),
      bp = vapply(idx, function(i) pgenlibr::GetVariantPos(pvar, i), integer(1)),
      a1 = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 2L), character(1)), # ALT, counted internally
      a2 = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 1L), character(1)), # REF
      z = seq_along(idx) / 10,
      stringsAsFactors = FALSE
    )
    input_file <- tempfile(fileext = ".txt")
    write.table(sumstats, input_file, quote = FALSE, row.names = FALSE)

    pop_wgt_df <- data.frame(
      pop = c("ACB", "ASW", "YRI", "LWK", "GBR", "FIN", "CHS", "CEU"),
      wgt = rep(1 / 8, 8),
      stringsAsFactors = FALSE
    )

    counts <- psam_pop_counts(psam_file, pop_col = "Population")
    stopifnot(all(c("pop", "n") %in% names(counts)))
    stopifnot(all(pop_wgt_df$pop %in% counts$pop))

    keep_dir <- tempfile("keep-files-")
    keep_paths <- write_keep_files_from_psam(
      psam_file,
      out_dir = keep_dir,
      pop_col = "Population",
      pops = pop_wgt_df$pop,
      include_fid = FALSE
    )
    stopifnot(all(file.exists(keep_paths)))
    stopifnot(all(names(keep_paths) %in% pop_wgt_df$pop))
    first_keep <- read.table(keep_paths[[1]], stringsAsFactors = FALSE)
    stopifnot(ncol(first_keep) == 1L)
    stopifnot(nrow(first_keep) > 0L)

    psam_tab <- read.table(psam_file, header = TRUE, comment.char = "", check.names = FALSE)
    names(psam_tab) <- sub("^#", "", names(psam_tab))

    # Default mode uses PSAM pedigree columns PAT/MAT and is ancestry-agnostic.
    ped_unrelated <- select_unrelated_from_psam(psam_file)
    stopifnot(!anyDuplicated(ped_unrelated$related_group))
    stopifnot(nrow(ped_unrelated) <= nrow(psam_tab))

    # Explicit relatedness clusters can also be supplied externally.
    related_meta <- data.frame(
      IID = psam_tab$IID,
      RELGROUP = paste0("RG", ceiling(seq_len(nrow(psam_tab)) / 2)),
      stringsAsFactors = FALSE
    )
    unrelated <- select_unrelated_from_psam(
      psam_file,
      group_col = "RELGROUP",
      sample_meta = related_meta,
      missing_group = "drop"
    )
    stopifnot(nrow(unrelated) == length(unique(related_meta$RELGROUP)))
    stopifnot(!anyDuplicated(unrelated$related_group))

    unrelated_keep <- file.path(tempfile("unrelated-keep-"), "unrelated.keep")
    written_unrelated <- write_unrelated_keep_file_from_psam(
      psam_file,
      keep_file = unrelated_keep,
      group_col = "RELGROUP",
      sample_meta = related_meta,
      include_fid = FALSE,
      missing_group = "drop"
    )
    stopifnot(file.exists(written_unrelated))
    unrelated_keep_tab <- read.table(written_unrelated, stringsAsFactors = FALSE)
    stopifnot(nrow(unrelated_keep_tab) == nrow(unrelated))
    stopifnot(ncol(unrelated_keep_tab) == 1L)

    ld <- computeLD_pgen(
      chr = 22,
      start_bp = min(sumstats$bp),
      end_bp = max(sumstats$bp),
      pop_wgt_df = pop_wgt_df,
      input_file = input_file,
      pgen_file = pgen_file,
      pvar_file = pvar_file,
      psam_file = psam_file,
      pop_col = "Population",
      af1_cutoff = -1
    )

    stopifnot(is.list(ld))
    stopifnot(all(c("snplist", "cormat") %in% names(ld)))
    stopifnot(is.data.frame(ld$snplist))
    stopifnot(all(c("rsid", "chr", "bp", "a1", "a2", "af1mix") %in% names(ld$snplist)))
    stopifnot(is.matrix(ld$cormat))
    stopifnot(nrow(ld$cormat) == nrow(ld$snplist))
    stopifnot(ncol(ld$cormat) == nrow(ld$snplist))
    stopifnot(max(abs(ld$cormat - t(ld$cormat)), na.rm = TRUE) < 1e-10)
    stopifnot(all(abs(diag(ld$cormat) - 1) < 1e-12))
    stopifnot(nrow(ld$snplist) > 10)
    stopifnot(all(ld$snplist$a1 %in% sumstats$a1)) # returned a1 is ALT

    measured_input <- sumstats[1:30, ]
    input_file2 <- tempfile(fileext = ".txt")
    write.table(measured_input, input_file2, quote = FALSE, row.names = FALSE)

    imp <- distmix_pgen(
      chr = 22,
      start_bp = pgenlibr::GetVariantPos(pvar, 1L),
      end_bp = pgenlibr::GetVariantPos(pvar, min(80L, nvar)),
      wing_size = 0,
      pop_wgt_df = pop_wgt_df,
      input_file = input_file2,
      pgen_file = pgen_file,
      pvar_file = pvar_file,
      psam_file = psam_file,
      pop_col = "Population",
      af1_cutoff = -1
    )

    stopifnot(is.data.frame(imp))
    stopifnot(all(c("rsid", "chr", "bp", "a1", "a2", "af1mix", "z", "pval", "info", "type") %in% names(imp)))
    stopifnot(nrow(imp) > 10)
    stopifnot(all(c(0L, 1L) %in% imp$type))
    stopifnot(all(is.finite(imp$z)))
    stopifnot(all(is.finite(imp$pval)))
    stopifnot(all(imp$pval >= 0 & imp$pval <= 1))
    stopifnot(all(is.finite(imp$info)))
    stopifnot(all(imp$info >= 0))

    # If input summary-stat a1/a2 is REF/ALT, GAUSS.pgen standardizes to
    # ALT/REF and flips the measured Z-score.
    flip_input <- sumstats
    flip_input$a1[1] <- sumstats$a2[1]
    flip_input$a2[1] <- sumstats$a1[1]
    flip_input$z[1] <- 9.25
    flip_file <- tempfile(fileext = ".txt")
    write.table(flip_input, flip_file, quote = FALSE, row.names = FALSE)
    matched <- GAUSS.pgen:::.gauss_match_sumstats_to_pvar(
      data.frame(
        chr = as.character(sumstats$chr[1]),
        bp = sumstats$bp[1],
        ref = sumstats$a2[1],
        alt = sumstats$a1[1],
        stringsAsFactors = FALSE
      ),
      GAUSS.pgen:::.gauss_read_sumstats_z(flip_file, 22, sumstats$bp[1], sumstats$bp[1])
    )
    stopifnot(identical(as.numeric(matched), -9.25))
  }
}
