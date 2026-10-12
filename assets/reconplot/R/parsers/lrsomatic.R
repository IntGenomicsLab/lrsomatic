## ---------------------------------------------------------------------------
## parsers/lrsomatic.R -- legacy ASCAT+Severus lrsomatic convenience wrapper.
##
## lrsomatic can run different CN/SV callers. Prefer naming those callers
## explicitly with --cn-source/--sv-source. This wrapper remains as the original
## ASCAT+Severus shortcut:
##
##   <sample_dir>/ascat/                        -> parsers/ascat.R    (CN)
##   <sample_dir>/variants/severus/somatic_SVs/ -> parsers/severus.R  (SVs)
##
## Equivalent to:
##   --cn-source ascat --cn-input <sample_dir>/ascat \
##   --sv-source severus --sv-input <sample_dir>/variants/severus
##
## For Wakhan CN from the same lrsomatic run, use:
##   --cn-source wakhan --sv-source severus --input <sample_dir>
##
## --input accepts the per-sample directory (output/P215003155), the pipeline
## outdir (output/, when it holds exactly one sample), or the run directory
## containing output/.
## ---------------------------------------------------------------------------

#' Resolve --input to an lrsomatic per-sample directory.
lrsomatic_sample_dir <- function(dir, sample = NULL) {
  if (is.null(dir)) stop("ASCAT+Severus lrsomatic shortcut needs --input")
  is_sample_dir <- function(d) dir.exists(file.path(d, "ascat")) ||
    dir.exists(file.path(d, "variants", "severus"))

  if (is_sample_dir(dir)) return(dir)

  ## descend through a run directory into output/
  for (d in c(file.path(dir, "output"), dir)) {
    if (!dir.exists(d)) next
    subs <- list.dirs(d, recursive = FALSE)
    subs <- subs[vapply(subs, is_sample_dir, logical(1))]
    if (!is.null(sample) && nzchar(sample)) {
      pinned <- subs[basename(subs) == sample | startsWith(basename(subs), sample)]
      if (length(pinned)) subs <- pinned
    }
    if (length(subs) == 1) return(subs[1])
    if (length(subs) > 1) {
      stop("Several lrsomatic sample directories under ", d, ": ",
           paste(basename(subs), collapse = ", "), ". Pass --sample or point ",
           "--input at one of them.")
    }
  }
  stop("Could not find an lrsomatic sample directory under ", dir,
       " (expected an 'ascat' and/or 'variants/severus' subdirectory)")
}

parse_lrsomatic <- function(args) {
  sample_dir <- lrsomatic_sample_dir(args$input, args$sample)
  log_msg("  lrsomatic sample directory: ", sample_dir)
  if (is.null(args$sample) || !nzchar(args$sample)) {
    args$sample <- basename(sample_dir)
    log_msg("  inferred sample: ", args$sample)
  }

  args_cn <- args; args_cn$input <- file.path(sample_dir, "ascat")
  args_sv <- args; args_sv$input <- file.path(sample_dir, "variants", "severus")
  if (!dir.exists(args_cn$input)) stop("No ascat/ directory under ", sample_dir)
  if (!dir.exists(args_sv$input)) stop("No variants/severus/ directory under ", sample_dir)

  ## ASCAT names its files "<sample>.tumour...", so do not pin on the bare
  ## sample id here; let the ascat parser infer its own prefix.
  args_cn$sample <- NULL
  p_cn <- parse_ascat(args_cn)
  p_sv <- parse_severus(args_sv)

  meta <- utils::modifyList(p_sv$meta %||% list(), p_cn$meta %||% list())
  meta$sample <- args$sample
  list(cn = p_cn$cn, sv = p_sv$sv, meta = meta, annotation_fn = p_cn$annotation_fn)
}

register_parser("lrsomatic", parse_lrsomatic,
                "Legacy lrsomatic shortcut: ASCAT copy number + Severus somatic SVs",
                provides = c("cn", "sv"))
