## ---------------------------------------------------------------------------
## plotting.R -- thin, caller-agnostic wrapper around ReConPlot::ReConPlot().
##
## Everything above this layer produces the two standard data frames; this file
## only decides sizing, titles, file names and per-panel error handling.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ReConPlot)
})

#' Default canvas size: grows with the number of side-by-side chromosome panels.
auto_plot_size <- function(regions, has_annotation = FALSE) {
  n <- nrow(regions)
  list(width = min(4 + 2.0 * n, 24),
       height = if (has_annotation) 4.6 else 3.4)
}

build_title <- function(user_title, meta, regions) {
  if (!is.null(user_title) && nzchar(user_title)) return(user_title)
  bits <- c()
  if (!is.null(meta$sample) && !is.na(meta$sample)) bits <- c(bits, meta$sample)
  if (!is.null(meta$purity) && !is.null(meta$ploidy)) {
    bits <- c(bits, sprintf("purity %.2f / ploidy %.2f",
                            as.numeric(meta$purity), as.numeric(meta$ploidy)))
  }
  chrs <- unique(regions$chr)
  where <- if (nrow(regions) >= 20 && all(regions$full)) {
    "genome-wide"
  } else if (length(chrs) > 6) {
    paste0(paste(chrs[1:6], collapse = ", "), ", +", length(chrs) - 6, " more")
  } else {
    paste(chrs, collapse = ", ")
  }
  bits <- c(bits, where)
  paste(bits, collapse = "  |  ")
}

#' Render one ReConPlot panel set.
#'
#' @param cn,sv validated data frames
#' @param regions chr_selection data frame
#' @param annotation optional data frame(chr, pos, y)
#' @param extra named list of further arguments forwarded to ReConPlot()
render_reconplot <- function(cn, sv, regions, title = "", genes = NULL,
                             annotation = NULL, extra = list()) {
  ## ReConPlot only reads chr/start/end; drop our bookkeeping columns.
  regions <- as.data.frame(regions[, c("chr", "start", "end")], stringsAsFactors = FALSE)
  rownames(regions) <- NULL
  call_args <- list(sv = sv, cnv = cn, chr_selection = regions, title = title)
  if (!is.null(genes) && length(genes)) call_args$genes <- genes
  if (!is.null(annotation) && nrow(annotation) > 0) {
    call_args$custom_annotation <- annotation
    call_args$ann_y_title <- "BAF"
    call_args$ann_one_scale <- TRUE
    call_args$ann_dot_size <- 0.15   # dense het-SNP track; keep dots small
  }
  call_args <- utils::modifyList(call_args, extra)
  do.call(ReConPlot::ReConPlot, call_args)
}

save_reconplot <- function(p, outdir, stem, width, height, formats = c("pdf")) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(0)
  for (fmt in formats) {
    path <- file.path(outdir, paste0(stem, ".", fmt))
    if (fmt == "png") {
      ggplot2::ggsave(path, plot = p, width = width, height = height,
                      units = "in", dpi = 300, limitsize = FALSE)
    } else {
      ggplot2::ggsave(path, plot = p, width = width, height = height,
                      units = "in", device = fmt, limitsize = FALSE)
    }
    paths <- c(paths, path)
  }
  paths
}

#' Render and save one plot per element of `region_sets`.
#'
#' @param region_sets list of chr_selection data frames
#' @return character vector of files written
plot_region_sets <- function(region_sets, cn, sv, meta, outdir, prefix,
                             title = NULL, genes = NULL, annotation_fn = NULL,
                             formats = c("pdf"), width = NULL, height = NULL,
                             extra = list()) {
  written <- character(0)
  for (i in seq_along(region_sets)) {
    regions <- region_sets[[i]]
    label <- region_label(regions)
    stem <- paste(c(prefix, label), collapse = "_")
    log_msg(sprintf("plotting %s (%d panel%s)", label, nrow(regions),
                    if (nrow(regions) == 1) "" else "s"))

    annotation <- if (!is.null(annotation_fn)) annotation_fn(regions) else NULL
    size <- auto_plot_size(regions, has_annotation = !is.null(annotation))
    w <- width %||% size$width
    h <- height %||% size$height

    p <- tryCatch(
      render_reconplot(cn, sv, regions,
                       title = build_title(title, meta, regions),
                       genes = genes, annotation = annotation, extra = extra),
      error = function(e) { log_msg("  ERROR: ", conditionMessage(e)); NULL })
    if (is.null(p)) next

    out <- tryCatch(save_reconplot(p, outdir, stem, w, h, formats),
                    error = function(e) { log_msg("  ERROR saving: ", conditionMessage(e)); character(0) })
    for (f in out) log_msg("  wrote ", f)
    written <- c(written, out)
  }
  written
}
