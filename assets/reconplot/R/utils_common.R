## ---------------------------------------------------------------------------
## utils_common.R -- caller-agnostic helpers shared by every ReConPlot parser.
##
## Nothing in this file knows about SAVANA (or any other caller). It defines the
## two data contracts that ReConPlot::ReConPlot() expects, plus the small amount
## of genome bookkeeping needed to turn a region string into a chr_selection
## data frame.
##
##   CN contract : chr, start, end, copyNumber, minorAlleleCopyNumber
##   SV contract : chr1, pos1, chr2, pos2, strands
##                 (strands is one of "+-", "-+", "++", "--", "TRA", "INS", "SBE")
## ---------------------------------------------------------------------------

MAIN_CHROMS <- paste0("chr", c(1:22, "X", "Y"))

## strand vocabulary understood by ReConPlot's colour mapper. Anything outside
## this set makes the package fail with an uninformative error, so we normalise
## aggressively on the way in.
VALID_STRANDS <- c("+-", "-+", "++", "--", "TRA", "INS", "SBE",
                   "DEL", "DUP", "h2hINV", "t2tINV")

CHROM_LENGTHS <- list(
  hg38 = c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
           chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
           chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
           chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
           chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
           chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415),
  hg19 = c(chr1 = 249250621, chr2 = 243199373, chr3 = 198022430, chr4 = 191154276,
           chr5 = 180915260, chr6 = 171115067, chr7 = 159138663, chr8 = 146364022,
           chr9 = 141213431, chr10 = 135534747, chr11 = 135006516, chr12 = 133851895,
           chr13 = 115169878, chr14 = 107349540, chr15 = 102531392, chr16 = 90354753,
           chr17 = 81195210, chr18 = 78077248, chr19 = 59128983, chr20 = 63025520,
           chr21 = 48129895, chr22 = 51304566, chrX = 155270560, chrY = 59373566)
)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), paste0(..., collapse = ""))
}

## --- chromosome naming ------------------------------------------------------

#' Coerce assorted chromosome spellings to the UCSC "chrN" style ReConPlot wants.
normalize_chrom <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(chr)?", "chr", x, ignore.case = TRUE)
  x <- sub("^chr(chr)+", "chr", x)
  x <- sub("^chrMT$", "chrM", x)
  x <- sub("^chr23$", "chrX", x)
  x <- sub("^chr24$", "chrY", x)
  x
}

is_main_chrom <- function(x) x %in% MAIN_CHROMS

#' Drop rows touching alt/random/decoy contigs, which ReConPlot rejects outright.
#'
#' @param df data frame
#' @param cols chromosome columns that must all be primary contigs
#' @param what label used in the log line
drop_nonstandard_chroms <- function(df, cols, what = "rows") {
  if (nrow(df) == 0) return(df)
  keep <- Reduce(`&`, lapply(cols, function(cc) is_main_chrom(df[[cc]])))
  n_drop <- sum(!keep)
  if (n_drop > 0) {
    bad <- unique(unlist(lapply(cols, function(cc) df[[cc]][!keep])))
    bad <- setdiff(bad, MAIN_CHROMS)
    log_msg(sprintf("  dropped %d %s on non-primary contigs (%s%s)",
                    n_drop, what, paste(utils::head(bad, 5), collapse = ", "),
                    if (length(bad) > 5) ", ..." else ""))
  }
  df[keep, , drop = FALSE]
}

## --- contract validation ----------------------------------------------------

validate_cn <- function(cn, drop_na_minor = FALSE) {
  req <- c("chr", "start", "end", "copyNumber", "minorAlleleCopyNumber")
  missing <- setdiff(req, names(cn))
  if (length(missing)) {
    stop("CN table is missing required column(s): ", paste(missing, collapse = ", "))
  }
  cn <- as.data.frame(cn, stringsAsFactors = FALSE)
  cn$chr <- normalize_chrom(cn$chr)
  cn$start <- as.numeric(cn$start)
  cn$end <- as.numeric(cn$end)
  cn$copyNumber <- suppressWarnings(as.numeric(cn$copyNumber))
  cn$minorAlleleCopyNumber <- suppressWarnings(as.numeric(cn$minorAlleleCopyNumber))

  cn <- drop_nonstandard_chroms(cn, "chr", "CN segments")
  bad <- is.na(cn$start) | is.na(cn$end) | is.na(cn$copyNumber)
  if (any(bad)) {
    log_msg(sprintf("  dropped %d CN segments with missing coordinates/total CN", sum(bad)))
    cn <- cn[!bad, , drop = FALSE]
  }
  if (drop_na_minor) {
    bad <- is.na(cn$minorAlleleCopyNumber)
    if (any(bad)) {
      log_msg(sprintf("  dropped %d CN segments with missing minor allele CN", sum(bad)))
      cn <- cn[!bad, , drop = FALSE]
    }
  }
  ## ReConPlot indexes cnv[, c("chr","start","end","copyNumber","minorAlleleCopyNumber")]
  ## positionally in places, so hand it exactly those columns in that order.
  cn <- cn[order(cn$chr, cn$start), req, drop = FALSE]
  rownames(cn) <- NULL
  if (nrow(cn) == 0) stop("No usable copy number segments after filtering.")
  cn
}

validate_sv <- function(sv, interchrom_as_tra = TRUE) {
  req <- c("chr1", "pos1", "chr2", "pos2", "strands")
  missing <- setdiff(req, names(sv))
  if (length(missing)) {
    stop("SV table is missing required column(s): ", paste(missing, collapse = ", "))
  }
  sv <- as.data.frame(sv, stringsAsFactors = FALSE)
  if (nrow(sv) == 0) return(sv[, union(req, names(sv)), drop = FALSE])

  sv$chr1 <- normalize_chrom(sv$chr1)
  sv$chr2 <- normalize_chrom(sv$chr2)
  sv$pos1 <- suppressWarnings(as.integer(round(as.numeric(sv$pos1))))
  sv$pos2 <- suppressWarnings(as.integer(round(as.numeric(sv$pos2))))
  sv$strands <- trimws(as.character(sv$strands))

  ## single breakends carry no mate; park them on their own locus so the
  ## chromosome filter below does not throw them away.
  sbe <- sv$strands == "SBE"
  if (any(sbe)) {
    sv$chr2[sbe] <- sv$chr1[sbe]
    sv$pos2[sbe] <- sv$pos1[sbe]
  }

  sv <- drop_nonstandard_chroms(sv, c("chr1", "chr2"), "SV junctions")
  bad <- is.na(sv$pos1) | is.na(sv$pos2) | is.na(sv$strands) | !(sv$strands %in% VALID_STRANDS)
  if (any(bad)) {
    log_msg(sprintf("  dropped %d SV junctions with unusable coordinates/orientation", sum(bad)))
    sv <- sv[!bad, , drop = FALSE]
  }
  if (interchrom_as_tra && nrow(sv) > 0) {
    tra <- sv$chr1 != sv$chr2 & !(sv$strands %in% c("INS", "SBE"))
    if (any(tra)) {
      sv$strands[tra] <- "TRA"
      log_msg(sprintf("  relabelled %d inter-chromosomal junctions as TRA", sum(tra)))
    }
  }
  rownames(sv) <- NULL
  sv
}

## --- region / chr_selection handling ---------------------------------------

#' Parse a region specification into a ReConPlot chr_selection data frame.
#'
#' Accepted forms (comma, semicolon or whitespace separated):
#'   "all"                       every primary contig present in the CN table
#'   "chr8"                      whole chromosome
#'   "chr8:120000000-130000000"  explicit window (commas/underscores allowed)
#'
#' @param spec character scalar, or NULL/"all"
#' @param cn validated CN table, used to bound whole-chromosome requests
#' @param genome one of names(CHROM_LENGTHS); NA lengths fall back to CN extent
parse_regions <- function(spec, cn, genome = "hg38") {
  lens <- CHROM_LENGTHS[[genome]]
  chrom_end <- function(chr) {
    from_ref <- if (!is.null(lens) && chr %in% names(lens)) unname(lens[[chr]]) else NA_real_
    from_cn <- suppressWarnings(max(cn$end[cn$chr == chr], na.rm = TRUE))
    if (!is.finite(from_cn)) from_cn <- NA_real_
    if (is.na(from_ref)) from_cn else from_ref
  }

  if (is.null(spec) || length(spec) == 0 || identical(tolower(trimws(spec)), "all")) {
    chrs <- MAIN_CHROMS[MAIN_CHROMS %in% unique(cn$chr)]
    if (length(chrs) == 0) stop("No primary contigs found in the CN table.")
    return(data.frame(chr = chrs, start = 0,
                      end = vapply(chrs, chrom_end, numeric(1)),
                      full = TRUE, stringsAsFactors = FALSE))
  }

  tokens <- unlist(strsplit(spec, "[,;[:space:]]+"))
  tokens <- tokens[nzchar(tokens)]
  out <- lapply(tokens, function(tok) {
    parts <- strsplit(tok, ":", fixed = TRUE)[[1]]
    chr <- normalize_chrom(parts[1])
    if (length(parts) == 1) {
      return(data.frame(chr = chr, start = 0, end = chrom_end(chr),
                        full = TRUE, stringsAsFactors = FALSE))
    }
    rng <- gsub("[,_]", "", parts[2])
    se <- strsplit(rng, "-", fixed = TRUE)[[1]]
    if (length(se) != 2) stop("Cannot parse region '", tok, "'. Use chr:start-end.")
    data.frame(chr = chr, start = as.numeric(se[1]), end = as.numeric(se[2]),
               full = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

read_regions_file <- function(path) {
  df <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                          comment.char = "#")
  if (ncol(df) < 3) stop("Regions BED file needs at least 3 columns: chr, start, end")
  data.frame(chr = normalize_chrom(df[[1]]),
             start = as.numeric(df[[2]]),
             end = as.numeric(df[[3]]),
             full = FALSE, stringsAsFactors = FALSE)
}

#' Keep only regions that ReConPlot can actually draw (primary contig + CN data).
sanitize_regions <- function(regions, cn) {
  regions <- as.data.frame(regions, stringsAsFactors = FALSE)
  regions$chr <- normalize_chrom(regions$chr)
  if (is.null(regions$full)) regions$full <- FALSE

  bad <- !is_main_chrom(regions$chr)
  if (any(bad)) {
    log_msg("  skipping unsupported contigs in region list: ",
            paste(unique(regions$chr[bad]), collapse = ", "))
    regions <- regions[!bad, , drop = FALSE]
  }
  missing_cn <- !(regions$chr %in% unique(cn$chr))
  if (any(missing_cn)) {
    log_msg("  skipping regions without copy number data: ",
            paste(unique(regions$chr[missing_cn]), collapse = ", "))
    regions <- regions[!missing_cn, , drop = FALSE]
  }
  if (nrow(regions) == 0) stop("No plottable regions left after filtering.")
  regions$start[is.na(regions$start) | regions$start < 0] <- 0
  regions <- regions[order(match(regions$chr, MAIN_CHROMS), regions$start), , drop = FALSE]
  rownames(regions) <- NULL
  regions
}

#' Short, filesystem-safe label for a region set (used in output file names).
region_label <- function(regions) {
  ## a full-chromosome sweep gets a name rather than a 24-part file stem
  if (nrow(regions) >= 20 && all(regions$full)) return("genome_wide")
  parts <- vapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i, ]
    if (isTRUE(r$full)) r$chr
    else sprintf("%s_%.1fMb-%.1fMb", r$chr, r$start / 1e6, r$end / 1e6)
  }, character(1))
  if (length(parts) > 6) parts <- c(parts[1:6], sprintf("and%dmore", length(parts) - 6))
  gsub("[^A-Za-z0-9._-]", "_", paste(parts, collapse = "_"))
}

#' Collapse repeated chromosomes into one spanning window.
#'
#' ReConPlot facets on `factor(chr, levels = unique(chr))`, so a multi-panel
#' selection containing the same chromosome twice fails with
#' "factor level [n] is duplicated". Merging is the useful behaviour: two
#' windows on one chromosome become the interval that spans both.
collapse_duplicate_chroms <- function(regions) {
  if (!anyDuplicated(regions$chr)) return(regions)
  dup <- unique(regions$chr[duplicated(regions$chr)])
  log_msg("  merging repeated chromosome(s) into a single panel: ",
          paste(dup, collapse = ", "))
  parts <- lapply(split(regions, regions$chr), function(g) {
    data.frame(chr = g$chr[1], start = min(g$start), end = max(g$end),
               full = any(g$full), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  out <- out[order(match(out$chr, MAIN_CHROMS)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
