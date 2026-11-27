## ---------------------------------------------------------------------------
## parsers/generic.R -- fallback parser for any caller not yet given its own file.
##
## Reads an explicit --cn-file and --sv-file and maps common column spellings
## onto the ReConPlot contract. Use it to sanity-check a new caller's output
## before writing a dedicated parser.
##
##   CN : chr/chromosome/seqnames, start, end, copyNumber (or nMajor+nMinor,
##        total_cn, CN), minorAlleleCopyNumber (or nMinor, minor_cn)
##   SV : headerless BEDPE (>=6 cols, orientation in col 9/10 or a "|"-packed
##        name column), or a headed table with chr1/pos1/chr2/pos2/strands.
## ---------------------------------------------------------------------------

.syn <- function(df, target, candidates) {
  if (target %in% names(df)) return(target)
  hit <- candidates[tolower(candidates) %in% tolower(names(df))]
  if (length(hit) == 0) return(NA_character_)
  names(df)[match(tolower(hit[1]), tolower(names(df)))]
}

generic_read_cn <- function(file) {
  df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE,
                          na.strings = c("NA", "nan", "NaN", ""))
  names(df) <- sub("^#", "", names(df))

  cchr   <- .syn(df, "chr", c("chromosome", "chrom", "seqnames", "CHR", "Chromosome"))
  cstart <- .syn(df, "start", c("startpos", "chromStart", "Start", "begin"))
  cend   <- .syn(df, "end", c("endpos", "chromEnd", "End", "stop"))
  ctot   <- .syn(df, "copyNumber", c("total_cn", "totalCN", "cn", "CN", "ntot", "tcn"))
  cmin   <- .syn(df, "minorAlleleCopyNumber", c("minor_cn", "minorCN", "nMinor",
                                                "nMin", "minor", "mcn"))
  cmaj   <- .syn(df, "nMajor", c("major_cn", "majorCN", "nMaj", "major"))

  for (nm in c(chr = cchr, start = cstart, end = cend)) {
    if (is.na(nm)) stop("Cannot find chr/start/end columns in ", basename(file))
  }
  out <- data.frame(chr = df[[cchr]], start = df[[cstart]], end = df[[cend]],
                    stringsAsFactors = FALSE)
  if (!is.na(ctot)) {
    out$copyNumber <- df[[ctot]]
  } else if (!is.na(cmaj) && !is.na(cmin)) {
    out$copyNumber <- as.numeric(df[[cmaj]]) + as.numeric(df[[cmin]])
  } else {
    stop("Cannot find a total copy number column in ", basename(file))
  }
  out$minorAlleleCopyNumber <- if (!is.na(cmin)) df[[cmin]] else NA_real_
  out
}

generic_read_sv <- function(file) {
  first <- readLines(file, n = 1L, warn = FALSE)
  headed <- grepl("chr1|chrom1|CHROM|chr_1", first, ignore.case = TRUE) &&
    !grepl("^chr[0-9XYM]", first)

  if (headed) {
    df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE)
    names(df) <- sub("^#", "", names(df))
    c1 <- .syn(df, "chr1", c("chrom1", "chromosome1", "chrA", "CHROM_A"))
    p1 <- .syn(df, "pos1", c("start1", "position1", "posA", "START_A"))
    c2 <- .syn(df, "chr2", c("chrom2", "chromosome2", "chrB", "CHROM_B"))
    p2 <- .syn(df, "pos2", c("start2", "position2", "posB", "START_B"))
    st <- .syn(df, "strands", c("strand", "orientation", "bp_notation", "BP_NOTATION"))
    if (any(is.na(c(c1, p1, c2, p2)))) {
      stop("Cannot find chr1/pos1/chr2/pos2 columns in ", basename(file))
    }
    strands <- if (!is.na(st)) df[[st]] else {
      s1 <- .syn(df, "strand1", c("str1", "orientation1"))
      s2 <- .syn(df, "strand2", c("str2", "orientation2"))
      if (any(is.na(c(s1, s2)))) stop("No strand information in ", basename(file))
      paste0(df[[s1]], df[[s2]])
    }
    sv <- data.frame(chr1 = df[[c1]], pos1 = df[[p1]], chr2 = df[[c2]],
                     pos2 = df[[p2]], strands = strands, stringsAsFactors = FALSE)
  } else {
    df <- data.table::fread(file, sep = "\t", header = FALSE, data.table = FALSE)
    if (ncol(df) < 6) stop("BEDPE ", basename(file), " has fewer than 6 columns")
    strands <- if (ncol(df) >= 10) {
      paste0(df[[9]], df[[10]])
    } else if (ncol(df) >= 7 && any(grepl("\\|", df[[7]]))) {
      vapply(strsplit(as.character(df[[7]]), "|", fixed = TRUE),
             function(p) p[length(p)], character(1))
    } else {
      stop("BEDPE ", basename(file), " carries no orientation (need cols 9-10 or a name field)")
    }
    sv <- data.frame(chr1 = df[[1]], pos1 = bedpe_point(df[[2]], df[[3]]),
                     chr2 = df[[4]], pos2 = bedpe_point(df[[5]], df[[6]]),
                     strands = strands, stringsAsFactors = FALSE)
  }
  sv$strands <- gsub("[<>]", "", trimws(sv$strands))
  sv$strands[toupper(sv$strands) == "INS"] <- "INS"
  sv$strands[toupper(sv$strands) %in% c("SBE", "SBND")] <- "SBE"
  sv
}

parse_generic <- function(args) {
  if (is.null(args$cn_file) || is.null(args$sv_file)) {
    stop("--source generic requires both --cn-file and --sv-file")
  }
  log_msg("  CN file: ", basename(args$cn_file))
  log_msg("  SV file: ", basename(args$sv_file))
  list(cn = generic_read_cn(args$cn_file),
       sv = generic_read_sv(args$sv_file),
       meta = list(sample = args$sample %||% NA_character_),
       annotation_fn = NULL)
}

register_parser("generic", parse_generic,
                "Any caller: explicit --cn-file/--sv-file with auto column mapping")
