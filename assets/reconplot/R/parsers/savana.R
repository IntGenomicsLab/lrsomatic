## ---------------------------------------------------------------------------
## parsers/savana.R -- turn a SAVANA output directory into ReConPlot inputs.
##
## Files consumed (all optional except the first two):
##   <sample>_segmented_absolute_copy_number.tsv  -> CN track
##   <sample>.classified.somatic.bedpe            -> SV track (default)
##   <sample>.classified.somatic.vcf              -> SV track (--sv-format vcf)
##   <sample>_fitted_purity_ploidy.tsv            -> purity/ploidy for the title
##   <sample>_allele_counts_hetSNPs.bed           -> optional BAF annotation track
##
## SAVANA writes breakpoint orientation as BP_NOTATION ("+-", "-+", "++", "--",
## "<INS>"), which is the same convention ReConPlot uses, so orientations pass
## through unchanged; only "<INS>"/"<SBE>" need renaming.
## ---------------------------------------------------------------------------

SAVANA_FILE_PATTERNS <- list(
  cn        = "_segmented_absolute_copy_number\\.tsv$",
  bedpe     = "\\.classified\\.somatic\\.bedpe$",
  vcf       = "\\.classified\\.somatic\\.vcf(\\.gz)?$",
  purity    = "_fitted_purity_ploidy\\.tsv$",
  hetsnp    = "_allele_counts_hetSNPs\\.bed(\\.gz)?$"
)

#' Locate a single SAVANA output file inside `dir`, optionally pinned to a sample.
savana_find_file <- function(dir, key, sample = NULL, required = TRUE) {
  pat <- SAVANA_FILE_PATTERNS[[key]]
  hits <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!is.null(sample) && nzchar(sample)) {
    pinned <- hits[startsWith(basename(hits), sample)]
    if (length(pinned)) hits <- pinned
  }
  if (length(hits) == 0) {
    if (required) stop("No SAVANA '", key, "' file matching /", pat, "/ in ", dir)
    return(NULL)
  }
  if (length(hits) > 1) {
    log_msg("  multiple '", key, "' files found; using ", basename(hits[1]),
            " (pass --sample to disambiguate)")
  }
  hits[1]
}

#' Infer the sample prefix from the copy number file name.
savana_infer_sample <- function(dir) {
  f <- list.files(dir, pattern = SAVANA_FILE_PATTERNS$cn)
  if (length(f) == 0) return(NA_character_)
  sub("_segmented_absolute_copy_number\\.tsv$", "", f[1])
}

## --- copy number ------------------------------------------------------------

savana_read_cn <- function(file) {
  df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE,
                          na.strings = c("NA", "nan", "NaN", ""))
  names(df)[names(df) == "chromosome"] <- "chr"
  need <- c("chr", "start", "end", "copyNumber", "minorAlleleCopyNumber")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("SAVANA CN file ", basename(file), " lacks column(s): ",
         paste(missing, collapse = ", "))
  }
  df[, need, drop = FALSE]
}

## --- structural variants ----------------------------------------------------

#' Normalise a SAVANA BP_NOTATION value to ReConPlot's `strands` vocabulary.
savana_notation_to_strands <- function(x) {
  x <- toupper(gsub("[<>]", "", trimws(as.character(x))))
  out <- rep(NA_character_, length(x))
  out[x %in% c("+-", "-+", "++", "--")] <- x[x %in% c("+-", "-+", "++", "--")]
  out[x == "INS"] <- "INS"
  out[x %in% c("SBE", "SBND", "BND_SINGLE")] <- "SBE"
  out
}

#' Pull the fields SAVANA packs into the BEDPE name column.
#'
#' The name looks like "ID_39590|1608bp|TUMOUR_12|--": an event id, the SV
#' length, the tumour read support and the breakpoint notation.
savana_split_bedpe_name <- function(name) {
  parts <- strsplit(as.character(name), "|", fixed = TRUE)
  pick <- function(p, rx, default = NA_character_) {
    m <- grep(rx, p, value = TRUE)
    if (length(m)) m[1] else default
  }
  data.frame(
    sv_id   = vapply(parts, function(p) p[1] %||% NA_character_, character(1)),
    svlen   = vapply(parts, function(p) {
      v <- pick(p, "^[0-9]+bp$"); suppressWarnings(as.numeric(sub("bp$", "", v)))
    }, numeric(1)),
    support = vapply(parts, function(p) {
      v <- pick(p, "^TUMOUR_[0-9]+$"); suppressWarnings(as.numeric(sub("^TUMOUR_", "", v)))
    }, numeric(1)),
    notation = vapply(parts, function(p) p[length(p)] %||% NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

#' BEDPE coordinates: SAVANA writes start == end == POS, while a spec-compliant
#' BEDPE writes start == POS-1. Handle both without guessing globally.
bedpe_point <- function(start, end) {
  ifelse(end - start == 1, end, start)
}

savana_read_sv_bedpe <- function(file, min_support = 0) {
  df <- data.table::fread(file, sep = "\t", header = FALSE, data.table = FALSE)
  if (ncol(df) < 7) stop("SAVANA BEDPE ", basename(file), " has fewer than 7 columns")
  names(df)[1:7] <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name")
  info <- savana_split_bedpe_name(df$name)

  sv <- data.frame(
    sv_id   = info$sv_id,
    chr1    = df$chrom1,
    pos1    = bedpe_point(df$start1, df$end1),
    chr2    = df$chrom2,
    pos2    = bedpe_point(df$start2, df$end2),
    strands = savana_notation_to_strands(info$notation),
    svlen   = info$svlen,
    support = info$support,
    stringsAsFactors = FALSE
  )
  savana_filter_sv(sv, min_support)
}

savana_read_sv_vcf <- function(file, min_support = 0, pass_only = TRUE) {
  vcf <- vcf_read_records(file, pass_only = pass_only)
  if (is.null(vcf) || nrow(vcf) == 0) return(savana_empty_sv())

  notation <- vcf_info_get(vcf$info, "BP_NOTATION")
  svlen    <- suppressWarnings(as.numeric(vcf_info_get(vcf$info, "SVLEN")))
  support  <- suppressWarnings(as.numeric(vcf_info_get(vcf$info, "TUMOUR_READ_SUPPORT")))
  mate     <- vcf_alt_mate(vcf$alt)

  sv <- data.frame(
    sv_id   = sub("_[12]$", "", vcf$id),
    chr1    = vcf$chrom,
    pos1    = vcf$pos,
    chr2    = ifelse(is.na(mate$chr), vcf$chrom, mate$chr),
    pos2    = ifelse(is.na(mate$pos), vcf$pos, mate$pos),
    strands = savana_notation_to_strands(notation),
    svlen   = svlen,
    support = support,
    stringsAsFactors = FALSE
  )
  sv <- dedupe_breakend_pairs(sv)   # BND records come in mate pairs
  savana_filter_sv(sv, min_support)
}

savana_empty_sv <- function() {
  data.frame(sv_id = character(), chr1 = character(), pos1 = numeric(),
             chr2 = character(), pos2 = numeric(), strands = character(),
             svlen = numeric(), support = numeric(), stringsAsFactors = FALSE)
}

savana_filter_sv <- function(sv, min_support = 0) {
  n0 <- nrow(sv)
  if (min_support > 0) {
    keep <- is.na(sv$support) | sv$support >= min_support
    if (any(!keep)) {
      log_msg(sprintf("  dropped %d SVs with tumour read support < %g",
                      sum(!keep), min_support))
    }
    sv <- sv[keep, , drop = FALSE]
  }
  log_msg(sprintf("  %d SV junctions read (%d after support filter)", n0, nrow(sv)))
  rownames(sv) <- NULL
  sv
}

## --- purity / ploidy --------------------------------------------------------

savana_read_purity <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  df <- utils::read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  as.list(df[1, , drop = FALSE])
}

## --- optional BAF annotation track -----------------------------------------

#' Build a `custom_annotation` data frame (chr, pos, y) of het-SNP BAF values.
#'
#' The het-SNP file is large, so it is read once, restricted to the primary
#' contigs and thinned; the per-region subset happens in the returned closure.
savana_baf_annotation_fn <- function(file, max_points_per_region = 5000) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  cache <- NULL
  function(regions) {
    if (is.null(cache)) {
      log_msg("  reading het-SNP BAF track: ", basename(file))
      d <- data.table::fread(file, sep = "\t", header = FALSE, data.table = FALSE,
                             select = c(1, 2, 12),
                             col.names = c("chr", "pos", "y"))
      d$chr <- normalize_chrom(d$chr)
      d <- d[is_main_chrom(d$chr) & !is.na(d$y), , drop = FALSE]
      cache <<- d
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

## --- entry point ------------------------------------------------------------

parse_savana <- function(args) {
  dir <- args$input
  sample <- args$sample
  if (is.null(sample) || !nzchar(sample)) {
    sample <- if (!is.null(dir)) savana_infer_sample(dir) else NA_character_
    if (!is.na(sample)) log_msg("  inferred sample: ", sample)
  }

  cn_file <- args$cn_file %||% savana_find_file(dir, "cn", sample)
  log_msg("  CN file: ", basename(cn_file))
  cn <- savana_read_cn(cn_file)

  sv_format <- match.arg(args$sv_format %||% "bedpe", c("bedpe", "vcf"))
  sv_file <- args$sv_file %||% savana_find_file(dir, sv_format, sample)
  log_msg("  SV file: ", basename(sv_file), " (", sv_format, ")")
  sv <- if (sv_format == "bedpe") {
    savana_read_sv_bedpe(sv_file, min_support = args$min_support %||% 0)
  } else {
    savana_read_sv_vcf(sv_file, min_support = args$min_support %||% 0)
  }

  meta <- list(sample = sample, cn_file = cn_file, sv_file = sv_file)
  pp <- savana_read_purity(args$purity_file %||%
                             (if (!is.null(dir)) savana_find_file(dir, "purity", sample, required = FALSE)))
  if (!is.null(pp)) {
    meta$purity <- pp$purity
    meta$ploidy <- pp$ploidy
    log_msg(sprintf("  purity=%s ploidy=%s", pp$purity, pp$ploidy))
  }

  ann_fn <- NULL
  if (isTRUE(args$baf_track)) {
    baf_file <- args$annotation_file %||%
      (if (!is.null(dir)) savana_find_file(dir, "hetsnp", sample, required = FALSE))
    ann_fn <- savana_baf_annotation_fn(baf_file,
                                       max_points_per_region = args$baf_max_points %||% 5000)
    if (is.null(ann_fn)) log_msg("  --baf-track requested but no het-SNP file found; skipping")
  }

  list(cn = cn, sv = sv, meta = meta, annotation_fn = ann_fn)
}

register_parser("savana", parse_savana,
                "SAVANA (segmented_absolute_copy_number.tsv + classified.somatic.bedpe/vcf)")
