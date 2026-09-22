## ---------------------------------------------------------------------------
## parsers/ascat.R -- ASCAT allele-specific copy number (CN component).
##
## Used by lrsomatic's long-read ASCAT step, but works for any ASCAT run
## (nf-core/sarek included) since the file names come from ASCAT itself.
##
## Files consumed:
##   <sample>.segments.txt      -> CN track (default; ascat.output$segments)
##   <sample>.segments_raw.txt  -> CN track with --ascat-cn-mode raw
##   <sample>.purityploidy.txt  -> purity (aberrant cell fraction) and ploidy
##   <sample>.tumour_tumourBAF.txt -> optional per-SNP BAF annotation track
##
## segments.txt holds the final fitted integer allele-specific calls
## (nMajor/nMinor), which is what ReConPlot wants: copyNumber = nMajor + nMinor
## and minorAlleleCopyNumber = nMinor. cnvs.txt is the same table minus the
## sample column, so it is only used as a fallback. segments_raw.txt is the
## pre-fit ASPCF segmentation and additionally carries the unrounded nAraw /
## nBraw, exposed via --ascat-cn-mode raw / raw-fractional.
##
## ASCAT writes chromosomes without a "chr" prefix; normalize_chrom() fixes that.
## ---------------------------------------------------------------------------

ASCAT_FILE_PATTERNS <- list(
  segments     = "\\.segments\\.txt$",
  segments_raw = "\\.segments_raw\\.txt$",
  cnvs         = "\\.cnvs\\.txt$",
  purityploidy = "\\.purityploidy\\.txt$",
  metrics      = "\\.metrics\\.txt$",
  baf          = "_tumourBAF\\.txt$"      # not _tumourBAF_rawBAF.txt
)

ascat_find_file <- function(dir, key, sample = NULL, required = TRUE) {
  pat <- ASCAT_FILE_PATTERNS[[key]]
  hits <- list.files(dir, pattern = pat, full.names = TRUE)
  if (key == "segments") hits <- hits[!grepl("segments_raw", hits)]
  if (!is.null(sample) && nzchar(sample)) {
    pinned <- hits[startsWith(basename(hits), sample)]
    if (length(pinned)) hits <- pinned
  }
  if (length(hits) == 0) {
    if (required) stop("No ASCAT '", key, "' file matching /", pat, "/ in ", dir)
    return(NULL)
  }
  if (length(hits) > 1) {
    log_msg("  multiple ASCAT '", key, "' files; using ", basename(hits[1]),
            " (pass --sample to disambiguate)")
  }
  hits[1]
}

ascat_infer_sample <- function(dir) {
  f <- list.files(dir, pattern = ASCAT_FILE_PATTERNS$segments)
  f <- f[!grepl("segments_raw", f)]
  if (length(f) == 0) return(NA_character_)
  sub("\\.segments\\.txt$", "", f[1])
}

#' Read ASCAT segments into the ReConPlot CN contract.
#'
#' @param mode "segments" (fitted integer calls), "raw" (pre-fit segmentation,
#'   still using the rounded nMajor/nMinor) or "raw-fractional" (pre-fit
#'   segmentation using the unrounded nAraw/nBraw)
ascat_read_cn <- function(file, mode = "segments", sample = NULL) {
  df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE,
                          na.strings = c("NA", "nan", "NaN", ""))
  need <- c("chr", "startpos", "endpos", "nMajor", "nMinor")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("ASCAT file ", basename(file), " lacks column(s): ", paste(missing, collapse = ", "))
  }
  ## segments.txt keeps a sample column ("<prefix>.tumour"); cnvs.txt does not.
  if ("sample" %in% names(df) && !is.null(sample) && nzchar(sample)) {
    hit <- grepl(sample, df$sample, fixed = TRUE)
    if (any(hit)) {
      if (any(!hit)) log_msg(sprintf("  kept %d/%d segments for sample '%s'",
                                     sum(hit), nrow(df), sample))
      df <- df[hit, , drop = FALSE]
    }
  }

  if (mode == "raw-fractional") {
    if (!all(c("nAraw", "nBraw") %in% names(df))) {
      stop("--ascat-cn-mode raw-fractional needs nAraw/nBraw (segments_raw.txt)")
    }
    total <- as.numeric(df$nAraw) + as.numeric(df$nBraw)
    minor <- as.numeric(df$nBraw)
  } else {
    total <- as.numeric(df$nMajor) + as.numeric(df$nMinor)
    minor <- as.numeric(df$nMinor)
  }

  data.frame(chr = df$chr, start = df$startpos, end = df$endpos,
             copyNumber = total, minorAlleleCopyNumber = minor,
             stringsAsFactors = FALSE)
}

ascat_read_purity <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  df <- utils::read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  purity <- if ("AberrantCellFraction" %in% names(df)) df$AberrantCellFraction[1] else NA_real_
  ploidy <- if ("Ploidy" %in% names(df)) df$Ploidy[1] else NA_real_
  if (is.na(purity) && is.na(ploidy)) return(NULL)
  list(purity = purity, ploidy = ploidy)
}

#' BAF annotation track from ASCAT's per-SNP tumour BAF table.
#'
#' Layout is: <rowname> Chromosome Position <sample>. Values are already
#' restricted to germline het SNPs, so no extra filtering is needed.
ascat_baf_annotation_fn <- function(file, max_points_per_region = 5000) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  cache <- NULL
  function(regions) {
    if (is.null(cache)) {
      log_msg("  reading ASCAT BAF track: ", basename(file))
      d <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE)
      cols <- names(d)
      chr_col <- cols[tolower(cols) == "chromosome"][1]
      pos_col <- cols[tolower(cols) == "position"][1]
      val_col <- setdiff(cols, c(chr_col, pos_col, cols[1]))
      if (is.na(chr_col) || is.na(pos_col) || length(val_col) == 0) {
        log_msg("  BAF file has an unexpected layout; skipping annotation track")
        return(NULL)
      }
      d <- data.frame(chr = normalize_chrom(d[[chr_col]]),
                      pos = as.numeric(d[[pos_col]]),
                      y = suppressWarnings(as.numeric(d[[val_col[1]]])),
                      stringsAsFactors = FALSE)
      cache <<- d[is_main_chrom(d$chr) & !is.na(d$y), , drop = FALSE]
    }
    out <- do.call(rbind, lapply(seq_len(nrow(regions)), function(i) {
      r <- regions[i, ]
      sub <- cache[cache$chr == r$chr & cache$pos >= r$start & cache$pos <= r$end, , drop = FALSE]
      if (nrow(sub) > max_points_per_region) {
        sub <- sub[sort(sample.int(nrow(sub), max_points_per_region)), , drop = FALSE]
      }
      sub
    }))
    if (is.null(out) || nrow(out) == 0) return(NULL)
    rownames(out) <- NULL
    out
  }
}

parse_ascat <- function(args) {
  dir <- args$input
  sample <- args$sample
  if ((is.null(sample) || !nzchar(sample)) && !is.null(dir)) {
    sample <- ascat_infer_sample(dir)
    if (!is.na(sample)) log_msg("  inferred ASCAT sample: ", sample)
  }

  mode <- match.arg(args$ascat_cn_mode %||% "segments",
                    c("segments", "raw", "raw-fractional"))
  cn_file <- args$cn_file
  if (is.null(cn_file)) {
    key <- if (mode == "segments") "segments" else "segments_raw"
    cn_file <- ascat_find_file(dir, key, sample, required = FALSE)
    if (is.null(cn_file) && mode == "segments") {
      cn_file <- ascat_find_file(dir, "cnvs", sample, required = FALSE)
      if (!is.null(cn_file)) log_msg("  segments.txt absent; falling back to cnvs.txt")
    }
    if (is.null(cn_file)) stop("No ASCAT copy number file found in ", dir)
  }
  log_msg("  CN file: ", basename(cn_file), " (mode: ", mode, ")")
  cn <- ascat_read_cn(cn_file, mode = mode, sample = sample)

  meta <- list(sample = sample, cn_file = cn_file, cn_caller = "ASCAT")
  pp <- ascat_read_purity(args$purity_file %||%
                            (if (!is.null(dir)) ascat_find_file(dir, "purityploidy", sample, FALSE)))
  if (!is.null(pp)) {
    meta$purity <- pp$purity
    meta$ploidy <- pp$ploidy
    log_msg(sprintf("  purity=%s ploidy=%s", pp$purity, round(as.numeric(pp$ploidy), 3)))
  }

  ann_fn <- NULL
  if (isTRUE(args$baf_track)) {
    baf_file <- args$annotation_file %||%
      (if (!is.null(dir)) ascat_find_file(dir, "baf", sample, required = FALSE))
    ann_fn <- ascat_baf_annotation_fn(baf_file,
                                      max_points_per_region = args$baf_max_points %||% 5000)
    if (is.null(ann_fn)) log_msg("  --baf-track requested but no ASCAT BAF file found; skipping")
  }

  list(cn = cn, meta = meta, annotation_fn = ann_fn)
}

register_parser("ascat", parse_ascat,
                "ASCAT allele-specific CN (segments.txt / segments_raw.txt)",
                provides = "cn")
