## ---------------------------------------------------------------------------
## parsers/wakhan.R -- Wakhan allele-specific copy number (CN component).
##
## Files consumed:
##   solutions_ranks.tsv                                      -> best solution
##   <solution>/bed_output/*_copynumbers_segments_HP_1.bed -> haplotype 1 CN
##   <solution>/bed_output/*_copynumbers_segments_HP_2.bed -> haplotype 2 CN
##
## The two haplotype BED files cover the full segmentation. ReConPlot wants
## total CN plus minor-allele CN, so join HP1/HP2 by segment coordinates and
## calculate total = HP1 + HP2, minor = min(HP1, HP2).
## ---------------------------------------------------------------------------

WAKHAN_FILE_PATTERNS <- list(
  ranks = "^solutions_ranks\\.tsv$",
  hp1   = "_copynumbers_segments_HP_1\\.bed(\\.gz)?$",
  hp2   = "_copynumbers_segments_HP_2\\.bed(\\.gz)?$"
)

wakhan_is_dir <- function(dir) {
  if (is.null(dir) || !dir.exists(dir)) return(FALSE)
  file.exists(file.path(dir, "solutions_ranks.tsv")) ||
    dir.exists(file.path(dir, "solution_1", "bed_output"))
}

#' Resolve --input to a Wakhan output directory.
wakhan_dir <- function(dir, sample = NULL) {
  if (is.null(dir)) stop("--cn-source wakhan needs --input or --cn-input")

  candidates <- c(dir, file.path(dir, "wakhan"))
  for (d in c(file.path(dir, "output"), dir)) {
    if (!dir.exists(d)) next
    subs <- list.dirs(d, recursive = FALSE)
    if (!is.null(sample) && nzchar(sample)) {
      pinned <- subs[basename(subs) == sample | startsWith(basename(subs), sample)]
      if (length(pinned)) subs <- pinned
    }
    candidates <- c(candidates, file.path(subs, "wakhan"), subs)
  }

  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  hits <- candidates[vapply(candidates, wakhan_is_dir, logical(1))]
  if (length(hits) == 1) return(hits[1])
  if (length(hits) > 1) {
    stop("Several Wakhan output directories found: ", paste(hits, collapse = ", "),
         ". Pass --sample or point --cn-input at one of them.")
  }
  stop("Could not find a Wakhan output directory under ", dir,
       " (expected solutions_ranks.tsv and/or solution_1/bed_output)")
}

wakhan_read_rank <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE,
                          na.strings = c("NA", "nan", "NaN", ""))
  if (nrow(df) == 0) return(NULL)
  if ("solution_rank" %in% names(df)) {
    df <- df[order(suppressWarnings(as.numeric(df$solution_rank))), , drop = FALSE]
  }
  df[1, , drop = FALSE]
}

wakhan_solution_dir <- function(dir, rank = NULL) {
  ranks <- wakhan_read_rank(file.path(dir, "solutions_ranks.tsv"))
  if (is.null(rank) && !is.null(ranks) && "solution_rank" %in% names(ranks)) {
    rank <- suppressWarnings(as.integer(ranks$solution_rank[1]))
  }
  rank <- rank %||% 1L
  sol_candidates <- character()
  if (!is.null(ranks) && "repository_name" %in% names(ranks)) {
    sol_candidates <- c(sol_candidates, file.path(dir, ranks$repository_name[1]))
  }
  sol_candidates <- c(sol_candidates, file.path(dir, paste0("solution_", rank)))
  subs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
  sol_candidates <- c(sol_candidates, subs[dir.exists(file.path(subs, "bed_output"))])
  sol_candidates <- unique(normalizePath(sol_candidates, mustWork = FALSE))
  sol_dir <- sol_candidates[dir.exists(file.path(sol_candidates, "bed_output"))]
  if (length(sol_dir) == 0) {
    stop("No Wakhan solution bed_output directory under ", dir,
         " (checked ", paste(basename(sol_candidates), collapse = ", "), ")")
  }
  sol_dir <- sol_dir[1]
  list(dir = sol_dir, rank = rank, ranks = ranks)
}

wakhan_find_file <- function(dir, key, sample = NULL, required = TRUE) {
  pat <- WAKHAN_FILE_PATTERNS[[key]]
  hits <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!is.null(sample) && nzchar(sample)) {
    pinned <- hits[startsWith(basename(hits), sample)]
    if (length(pinned)) hits <- pinned
  }
  if (length(hits) == 0) {
    if (required) stop("No Wakhan '", key, "' file matching /", pat, "/ in ", dir)
    return(NULL)
  }
  if (length(hits) > 1) {
    log_msg("  multiple Wakhan '", key, "' files; using ", basename(hits[1]),
            " (pass --sample to disambiguate)")
  }
  hits[1]
}

wakhan_read_hp_bed <- function(file) {
  df <- data.table::fread(file, sep = "\t", header = TRUE, data.table = FALSE,
                          skip = "#chr\tstart\tend", na.strings = c("NA", "nan", "NaN", ""))
  names(df) <- sub("^#", "", names(df))
  need <- c("chr", "start", "end", "copynumber_state")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("Wakhan BED file ", basename(file), " lacks column(s): ",
         paste(missing, collapse = ", "))
  }
  data.frame(chr = df$chr, start = df$start, end = df$end,
             cn = suppressWarnings(as.numeric(df$copynumber_state)),
             stringsAsFactors = FALSE)
}

wakhan_read_cn <- function(hp1_file, hp2_file) {
  hp1 <- wakhan_read_hp_bed(hp1_file)
  hp2 <- wakhan_read_hp_bed(hp2_file)
  cn <- merge(hp1, hp2, by = c("chr", "start", "end"), suffixes = c("1", "2"),
              all = FALSE, sort = FALSE)
  if (nrow(cn) == 0) stop("Wakhan HP1/HP2 BED files share no segment coordinates")
  if (nrow(cn) < max(nrow(hp1), nrow(hp2))) {
    log_msg(sprintf("  kept %d shared Wakhan HP segments from %d/%d rows",
                    nrow(cn), nrow(hp1), nrow(hp2)))
  }
  data.frame(chr = cn$chr, start = cn$start, end = cn$end,
             copyNumber = cn$cn1 + cn$cn2,
             minorAlleleCopyNumber = pmin(cn$cn1, cn$cn2),
             stringsAsFactors = FALSE)
}

parse_wakhan <- function(args) {
  dir <- wakhan_dir(args$input, args$sample)
  sol <- wakhan_solution_dir(dir)
  bed_dir <- file.path(sol$dir, "bed_output")
  if (!dir.exists(bed_dir)) stop("No Wakhan bed_output directory under ", sol$dir)

  hp1 <- args$cn_file %||% wakhan_find_file(bed_dir, "hp1", args$sample)
  hp2 <- wakhan_find_file(bed_dir, "hp2", args$sample)
  log_msg("  Wakhan directory: ", dir)
  log_msg("  Wakhan solution: ", basename(sol$dir))
  log_msg("  CN files: ", basename(hp1), " + ", basename(hp2))
  cn <- wakhan_read_cn(hp1, hp2)

  meta <- list(sample = args$sample %||% NA_character_, cn_file = hp1, cn_caller = "Wakhan",
               wakhan_solution = basename(sol$dir))
  if (!is.null(sol$ranks)) {
    meta$sample <- meta$sample %||% sub("_.*$", "", basename(hp1))
    if ("cell_purity" %in% names(sol$ranks)) meta$purity <- sol$ranks$cell_purity[1]
    if ("ploidy" %in% names(sol$ranks)) meta$ploidy <- sol$ranks$ploidy[1]
    if ("confidence" %in% names(sol$ranks)) meta$wakhan_confidence <- sol$ranks$confidence[1]
    log_msg(sprintf("  purity=%s ploidy=%s", meta$purity %||% NA, meta$ploidy %||% NA))
  }

  list(cn = cn, meta = meta, annotation_fn = NULL)
}

register_parser("wakhan", parse_wakhan,
                "Wakhan allele-specific CN (<solution>/bed_output HP BEDs)",
                provides = "cn")
