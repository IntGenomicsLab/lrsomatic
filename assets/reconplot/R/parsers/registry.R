## ---------------------------------------------------------------------------
## registry.R -- plug-in table mapping a caller name to a parser function.
##
## Adding support for a new SV/CN caller means dropping one file into
## R/parsers/ that ends with a register_parser() call. Nothing else changes.
##
## A parser is a function(args) -> list(...) where `args` is the list of CLI
## options (see run_reconplot.R). What it must return depends on `provides`:
##
##   provides = c("cn", "sv")  full parser: list(cn = <df>, sv = <df>, ...)
##   provides = "cn"           CN component: list(cn = <df>, ...)
##   provides = "sv"           SV component: list(sv = <df>, ...)
##
## Optional extra elements: `meta` (named list; sample/purity/ploidy feed the
## plot title) and `annotation_fn` (function(regions) -> data.frame(chr,pos,y)).
##
## Components exist because a pipeline may call CN and SVs with different tools
## (lrsomatic: ASCAT + Severus). They compose via --cn-source / --sv-source, so
## any CN caller can be paired with any SV caller.
##
## Parsers should return raw-ish tables; validate_cn()/validate_sv() from
## utils_common.R are applied centrally afterwards.
## ---------------------------------------------------------------------------

.parser_registry <- new.env(parent = emptyenv())

register_parser <- function(name, fn, description = "", provides = c("cn", "sv")) {
  stopifnot(is.character(name), length(name) == 1, is.function(fn))
  provides <- match.arg(provides, c("cn", "sv"), several.ok = TRUE)
  assign(name, list(fn = fn, description = description, provides = provides),
         envir = .parser_registry)
  invisible(NULL)
}

get_parser_entry <- function(name) {
  if (is.null(name) || !nzchar(name) ||
      !exists(name, envir = .parser_registry, inherits = FALSE)) {
    stop("Unknown source '", name, "'. Available: ",
         paste(list_parsers()$name, collapse = ", "))
  }
  get(name, envir = .parser_registry, inherits = FALSE)
}

get_parser <- function(name) get_parser_entry(name)$fn

parser_provides <- function(name) get_parser_entry(name)$provides

list_parsers <- function() {
  nms <- sort(ls(.parser_registry))
  data.frame(
    name = nms,
    provides = vapply(nms, function(n)
      paste(get(n, envir = .parser_registry)$provides, collapse = "+"), character(1)),
    description = vapply(nms, function(n)
      get(n, envir = .parser_registry)$description, character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Source every parser implementation found in `dir`.
load_parsers <- function(dir) {
  files <- setdiff(list.files(dir, pattern = "\\.[Rr]$", full.names = TRUE),
                   file.path(dir, "registry.R"))
  for (f in sort(files)) sys.source(f, envir = globalenv())
  invisible(list_parsers())
}

#' Work out which parser supplies CN and which supplies SVs, then run them.
#'
#' Single-source mode calls one full parser. Split mode (--cn-source and/or
#' --sv-source) calls one parser for each half, each with its own --cn-input /
#' --sv-input if given, and merges the results.
run_parsers <- function(args) {
  cn_src <- args$cn_source %||% args$source
  sv_src <- args$sv_source %||% args$source
  split_mode <- !identical(cn_src, sv_src) ||
    !is.null(args$cn_input) || !is.null(args$sv_input)

  if (!split_mode) {
    prov <- parser_provides(cn_src)
    if (!all(c("cn", "sv") %in% prov)) {
      stop("Source '", cn_src, "' only provides '", paste(prov, collapse = "+"),
           "'. Pair it with another via --cn-source/--sv-source, e.g. ",
           "--cn-source ascat --sv-source severus")
    }
    log_msg("source: ", cn_src)
    return(get_parser(cn_src)(args))
  }

  log_msg("CN source: ", cn_src, "   SV source: ", sv_src)
  for (s in c(cn = cn_src, sv = sv_src)) invisible(get_parser_entry(s))

  args_cn <- args; args_cn$input <- args$cn_input %||% args$input
  args_sv <- args; args_sv$input <- args$sv_input %||% args$input

  p_cn <- get_parser(cn_src)(args_cn)
  p_sv <- get_parser(sv_src)(args_sv)
  if (is.null(p_cn$cn)) stop("Source '", cn_src, "' returned no copy number table")
  if (is.null(p_sv$sv)) stop("Source '", sv_src, "' returned no SV table")

  list(cn = p_cn$cn, sv = p_sv$sv,
       ## CN metadata wins: purity/ploidy come from the CN fit
       meta = utils::modifyList(p_sv$meta %||% list(), p_cn$meta %||% list()),
       annotation_fn = p_cn$annotation_fn %||% p_sv$annotation_fn)
}
