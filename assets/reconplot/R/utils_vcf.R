## ---------------------------------------------------------------------------
## utils_vcf.R -- minimal, caller-agnostic VCF reading for SV parsers.
##
## Deliberately not a full VCF library: SV callers only need the fixed columns,
## a few INFO keys and the breakend ALT notation. Shared by savana.R and
## severus.R so both agree on how a BND mate is decoded.
## ---------------------------------------------------------------------------

#' Read a (optionally bgzipped) VCF into a data frame of the fixed columns.
#'
#' @return data frame with chrom, pos, id, ref, alt, qual, filter, info,
#'         format, sample1 (the first sample column, if present)
vcf_read_records <- function(file, pass_only = TRUE) {
  con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "##")]
  if (length(lines) < 2) return(NULL)

  f <- strsplit(lines[-1], "\t", fixed = TRUE)
  ncol_max <- max(lengths(f))
  f <- do.call(rbind, lapply(f, function(x) c(x, rep(NA_character_, ncol_max - length(x)))))

  out <- data.frame(
    chrom  = f[, 1], pos = suppressWarnings(as.numeric(f[, 2])), id = f[, 3],
    ref    = f[, 4], alt = f[, 5], qual = f[, 6], filter = f[, 7], info = f[, 8],
    format = if (ncol_max >= 9) f[, 9] else NA_character_,
    sample1 = if (ncol_max >= 10) f[, 10] else NA_character_,
    stringsAsFactors = FALSE)

  if (pass_only) {
    keep <- out$filter %in% c("PASS", ".") | is.na(out$filter)
    if (any(!keep)) log_msg(sprintf("  dropped %d non-PASS VCF records", sum(!keep)))
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Extract one INFO key. Returns NA where the key is absent.
vcf_info_get <- function(info, key) {
  rx <- paste0("(^|;)", key, "=([^;]*)")
  m <- regexpr(rx, info)
  out <- rep(NA_character_, length(info))
  hit <- which(m > 0)
  if (length(hit)) {
    out[hit] <- sub(paste0("^.*?", key, "="), "", regmatches(info, m), perl = TRUE)
  }
  out
}

#' TRUE where a valueless INFO flag (e.g. PRECISE, IMPRECISE) is present.
vcf_info_flag <- function(info, key) {
  grepl(paste0("(^|;)", key, "(;|$)"), info)
}

#' Decode the mate locus from a breakend ALT allele: t[chr:pos[ or ]chr:pos]t
#'
#' @return list(chr, pos); both NA where the ALT is not a breakend
vcf_alt_mate <- function(alt) {
  rx <- "[][][^][]+:[0-9]+[][]"
  m <- regexpr(rx, alt)
  chr <- rep(NA_character_, length(alt)); pos <- rep(NA_real_, length(alt))
  hit <- which(m > 0)
  if (length(hit)) {
    clean <- gsub("[][]", "", regmatches(alt, m))
    chr[hit] <- sub(":[0-9]+$", "", clean)
    pos[hit] <- as.numeric(sub("^.*:", "", clean))
  }
  list(chr = chr, pos = pos)
}

#' Pull one FORMAT subfield (e.g. "DV") out of the sample column.
vcf_format_get <- function(format, sample, key) {
  vapply(seq_along(format), function(i) {
    if (is.na(format[i]) || is.na(sample[i])) return(NA_character_)
    keys <- strsplit(format[i], ":", fixed = TRUE)[[1]]
    j <- match(key, keys)
    if (is.na(j)) return(NA_character_)
    vals <- strsplit(sample[i], ":", fixed = TRUE)[[1]]
    if (j > length(vals)) NA_character_ else vals[j]
  }, character(1))
}

#' Collapse mate-paired breakend records to one row per junction.
#'
#' Keyed on the unordered breakpoint pair so it works whether or not the caller
#' emits MATE_ID.
dedupe_breakend_pairs <- function(sv) {
  if (nrow(sv) == 0) return(sv)
  a <- paste0(sv$chr1, ":", sv$pos1)
  b <- paste0(sv$chr2, ":", sv$pos2)
  key <- ifelse(a < b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))
  sv[!duplicated(key), , drop = FALSE]
}
