## ---------------------------------------------------------------------------
## parsers/severus.R -- Severus long-read somatic SV calls (SV component).
##
## Files consumed:
##   somatic_SVs/severus_somatic.vcf.gz         -> SV track (the full callset)
##   somatic_SVs/breakpoint_clusters_list.tsv   -> optional cluster summary log
##
## Why the VCF and not breakpoint_clusters.tsv: the cluster file lists only the
## junctions Severus assigned to a complex-rearrangement cluster. On sample
## P215003155 that is 213 of the 330 junctions in the VCF (every cluster
## junction is present in the VCF; 117 VCF junctions are in no cluster). It is
## a filter/annotation layer, not the callset. Cluster membership is available
## from the VCF anyway as INFO/CLUSTERID, so --cluster-id filters without ever
## reading the TSV.
##
## Severus INFO/STRANDS uses the same orientation convention as ReConPlot and
## SAVANA (DEL "+-", DUP "-+", INV "++"/"--"), verified against records both
## callers report for this sample, so orientations pass straight through.
##
## Note on defaults: Severus runs at --min-support 3 in lrsomatic and emits
## small VNTR indels alongside rearrangements (99 junctions under 1 kb here,
## versus 19 from SAVANA). Nothing is filtered by default -- use --min-svlen
## and --exclude-vntr to get a rearrangement-only view.
## ---------------------------------------------------------------------------

SEVERUS_FILE_PATTERNS <- list(
  somatic_vcf  = "^severus_somatic\\.vcf(\\.gz)?$",
  all_vcf      = "^severus_all\\.vcf(\\.gz)?$",
  clusters     = "^breakpoint_clusters\\.tsv$",
  cluster_list = "^breakpoint_clusters_list\\.tsv$"
)

#' Find a Severus file, searching the given directory and the usual subdirs.
#'
#' Accepts the severus/ directory, the somatic_SVs/ directory, or an lrsomatic
#' per-sample directory, so callers do not have to know the layout.
severus_find_file <- function(dir, key, required = TRUE, somatic = TRUE) {
  pat <- SEVERUS_FILE_PATTERNS[[key]]
  sub <- if (somatic) "somatic_SVs" else "all_SVs"
  candidates <- c(dir,
                  file.path(dir, sub),
                  file.path(dir, "severus", sub),
                  file.path(dir, "variants", "severus", sub))
  for (d in candidates) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = pat, full.names = TRUE)
    if (length(hits)) return(hits[1])
  }
  if (required) stop("No Severus '", key, "' file matching /", pat, "/ under ", dir)
  NULL
}

#' Log the largest Severus clusters so the user can pick one for --cluster-id.
severus_log_cluster_list <- function(file, top = 5) {
  if (is.null(file) || !file.exists(file)) return(invisible(NULL))
  df <- utils::read.table(file, sep = "\t", header = FALSE, skip = 1,
                          stringsAsFactors = FALSE, quote = "", fill = TRUE)
  if (ncol(df) < 5) return(invisible(NULL))
  names(df)[c(1, 2, 5)] <- c("cluster_id", "type", "sv_count")
  df <- df[order(-suppressWarnings(as.numeric(df$sv_count))), , drop = FALSE]
  n <- min(top, nrow(df))
  log_msg(sprintf("  largest Severus clusters: %s",
                  paste(sprintf("%s (%s, %s SVs)", df$cluster_id[1:n],
                                df$type[1:n], df$sv_count[1:n]), collapse = "; ")))
  invisible(df)
}

#' Read severus_somatic.vcf.gz into the ReConPlot SV contract.
severus_read_vcf <- function(file, min_support = 0, min_svlen = 0,
                             exclude_vntr = FALSE, cluster_id = NULL,
                             clustered_only = FALSE, pass_only = TRUE) {
  vcf <- vcf_read_records(file, pass_only = pass_only)
  if (is.null(vcf) || nrow(vcf) == 0) return(severus_empty_sv())

  svtype  <- vcf_info_get(vcf$info, "SVTYPE")
  strands <- vcf_info_get(vcf$info, "STRANDS")
  end     <- suppressWarnings(as.numeric(vcf_info_get(vcf$info, "END")))
  svlen   <- suppressWarnings(as.numeric(vcf_info_get(vcf$info, "SVLEN")))
  cluster <- vcf_info_get(vcf$info, "CLUSTERID")
  detail  <- vcf_info_get(vcf$info, "DETAILED_TYPE")
  vntr    <- !is.na(vcf_info_get(vcf$info, "INSIDE_VNTR"))
  support <- suppressWarnings(as.numeric(vcf_format_get(vcf$format, vcf$sample1, "DV")))

  ## Mate locus: breakends carry it in the ALT allele, symbolic ALTs
  ## (<DEL>/<DUP>/<INV>) carry it as INFO/END on the same chromosome.
  mate <- vcf_alt_mate(vcf$alt)
  chr2 <- ifelse(!is.na(mate$chr), mate$chr, vcf$chrom)
  pos2 <- ifelse(!is.na(mate$pos), mate$pos,
                 ifelse(!is.na(end), end, vcf$pos))

  ## STRANDS is absent on insertions; fall back to the SV type.
  strands[is.na(strands) & svtype == "INS"] <- "INS"
  strands[is.na(strands) & svtype == "DEL"] <- "+-"
  strands[is.na(strands) & svtype == "DUP"] <- "-+"

  sv <- data.frame(
    sv_id    = sub("_[12]$", "", vcf$id),
    chr1     = vcf$chrom, pos1 = vcf$pos,
    chr2     = chr2,      pos2 = pos2,
    strands  = strands,
    svtype   = svtype,
    detailed_type = detail,
    svlen    = svlen,
    support  = support,
    cluster_id = cluster,
    inside_vntr = vntr,
    stringsAsFactors = FALSE)

  n_records <- nrow(sv)
  sv <- dedupe_breakend_pairs(sv)
  log_msg(sprintf("  %d VCF records -> %d junctions after mate collapsing",
                  n_records, nrow(sv)))
  severus_filter_sv(sv, min_support = min_support, min_svlen = min_svlen,
                    exclude_vntr = exclude_vntr, cluster_id = cluster_id,
                    clustered_only = clustered_only)
}

severus_empty_sv <- function() {
  data.frame(sv_id = character(), chr1 = character(), pos1 = numeric(),
             chr2 = character(), pos2 = numeric(), strands = character(),
             svtype = character(), detailed_type = character(), svlen = numeric(),
             support = numeric(), cluster_id = character(), inside_vntr = logical(),
             stringsAsFactors = FALSE)
}

severus_filter_sv <- function(sv, min_support = 0, min_svlen = 0,
                              exclude_vntr = FALSE, cluster_id = NULL,
                              clustered_only = FALSE) {
  drop <- function(sv, keep, why) {
    if (any(!keep)) log_msg(sprintf("  dropped %d junctions: %s", sum(!keep), why))
    sv[keep, , drop = FALSE]
  }
  if (min_support > 0) {
    sv <- drop(sv, is.na(sv$support) | sv$support >= min_support,
               sprintf("read support < %g", min_support))
  }
  if (min_svlen > 0) {
    ## A length only means something within a chromosome; never use it to
    ## discard translocations, which carry no SVLEN.
    intra <- sv$chr1 == sv$chr2
    span <- ifelse(!is.na(sv$svlen), abs(sv$svlen), abs(sv$pos2 - sv$pos1))
    sv <- drop(sv, !intra | is.na(span) | span >= min_svlen,
               sprintf("intra-chromosomal span < %g bp", min_svlen))
  }
  if (exclude_vntr) {
    sv <- drop(sv, !sv$inside_vntr, "inside a VNTR")
  }
  if (!is.null(cluster_id) && nzchar(cluster_id)) {
    wanted <- trimws(unlist(strsplit(cluster_id, ",")))
    sv <- drop(sv, !is.na(sv$cluster_id) & sv$cluster_id %in% wanted,
               paste("not in cluster", paste(wanted, collapse = "/")))
  } else if (clustered_only) {
    sv <- drop(sv, !is.na(sv$cluster_id), "not assigned to any cluster")
  }
  rownames(sv) <- NULL
  sv
}

parse_severus <- function(args) {
  dir <- args$input
  somatic <- !isTRUE(args$severus_all)
  sv_file <- args$sv_file %||%
    severus_find_file(dir, if (somatic) "somatic_vcf" else "all_vcf", somatic = somatic)
  log_msg("  SV file: ", basename(sv_file),
          if (somatic) " (somatic)" else " (all SVs, germline included)")

  sv <- severus_read_vcf(
    sv_file,
    min_support    = args$min_support %||% 0,
    min_svlen      = args$min_svlen %||% 0,
    exclude_vntr   = isTRUE(args$exclude_vntr),
    cluster_id     = args$cluster_id,
    clustered_only = isTRUE(args$clustered_only))

  if (nrow(sv) > 0) {
    tab <- table(sv$svtype, useNA = "no")
    log_msg("  Severus SV types: ",
            paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = "  "))
  }
  severus_log_cluster_list(severus_find_file(dir, "cluster_list", required = FALSE,
                                             somatic = somatic))

  list(sv = sv, meta = list(sample = args$sample, sv_file = sv_file,
                            sv_caller = "Severus"))
}

register_parser("severus", parse_severus,
                "Severus long-read SVs (somatic_SVs/severus_somatic.vcf.gz)",
                provides = "sv")
