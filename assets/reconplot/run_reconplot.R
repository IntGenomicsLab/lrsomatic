#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## run_reconplot.R -- command line front end for ReConPlot.
##
## Pipeline:  <caller output> --[parser]--> {cn, sv} --[validate]--> ReConPlot
##
## The only caller-specific code lives in R/parsers/. Everything here is
## generic, so supporting a new caller means adding one parser file.
##
## Examples
##   ./run_reconplot.R --cn-source savana --sv-source savana \
##       --input  /path/to/savana/sample \
##       --outdir /path/to/ReConPlot_output
##
##   ./run_reconplot.R --cn-source wakhan --sv-source severus --input DIR --outdir OUT \
##       --regions "chr8,chr17:0-30000000" --layout together --genes MYC,TP53
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

## --- locate our own installation directory ---------------------------------
script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) return(normalizePath(dirname(f[1])))
  normalizePath(".")
}
SCRIPT_DIR <- script_dir()
source(file.path(SCRIPT_DIR, "R", "utils_common.R"))
source(file.path(SCRIPT_DIR, "R", "utils_vcf.R"))
source(file.path(SCRIPT_DIR, "R", "parsers", "registry.R"))
source(file.path(SCRIPT_DIR, "R", "plotting.R"))
load_parsers(file.path(SCRIPT_DIR, "R", "parsers"))

## --- options ----------------------------------------------------------------
option_list <- list(
  make_option("--source", type = "character", default = "savana",
              help = "Single parser for CN+SVs; prefer --cn-source/--sv-source for mixed callers [default %default]"),
  make_option("--cn-source", type = "character", default = NULL,
              help = "Parser for copy number only (e.g. ascat); overrides --source"),
  make_option("--sv-source", type = "character", default = NULL,
              help = "Parser for SVs only (e.g. severus); overrides --source"),
  make_option("--cn-input", type = "character", default = NULL,
              help = "Directory for the CN source [default: --input]"),
  make_option("--sv-input", type = "character", default = NULL,
              help = "Directory for the SV source [default: --input]"),
  make_option("--list-sources", action = "store_true", default = FALSE,
              help = "Print the registered parsers and exit"),
  make_option("--input", type = "character", default = NULL,
              help = "Caller output directory (parser discovers files inside)"),
  make_option("--sample", type = "character", default = NULL,
              help = "Sample prefix; inferred from file names when omitted"),
  make_option("--cn-file", type = "character", default = NULL,
              help = "Explicit copy number file (overrides discovery)"),
  make_option("--sv-file", type = "character", default = NULL,
              help = "Explicit SV file (overrides discovery)"),
  make_option("--purity-file", type = "character", default = NULL,
              help = "Explicit purity/ploidy file"),
  make_option("--sv-format", type = "character", default = "bedpe",
              help = "SAVANA SV input: bedpe or vcf [default %default]"),
  make_option("--min-support", type = "double", default = 0,
              help = "Drop SVs with tumour read support below this [default %default]"),
  make_option("--min-svlen", type = "double", default = 0,
              help = paste("Drop intra-chromosomal SVs shorter than this many bp;",
                           "translocations are never dropped [default %default]")),
  make_option("--exclude-vntr", action = "store_true", default = FALSE,
              help = "Drop SVs flagged as inside a VNTR (Severus)"),
  make_option("--cluster-id", type = "character", default = NULL,
              help = "Keep only SVs in these Severus cluster IDs (comma-separated)"),
  make_option("--clustered-only", action = "store_true", default = FALSE,
              help = "Keep only SVs assigned to some Severus cluster"),
  make_option("--severus-all", action = "store_true", default = FALSE,
              help = "Use severus_all.vcf.gz (germline included) instead of the somatic set"),
  make_option("--ascat-cn-mode", type = "character", default = "segments",
              help = paste("ASCAT CN table: segments (fitted integer calls),",
                           "raw (pre-fit segmentation), raw-fractional",
                           "(pre-fit, unrounded nAraw/nBraw) [default %default]")),

  make_option("--outdir", type = "character", default = NULL,
              help = "Output directory (required)"),
  make_option("--prefix", type = "character", default = NULL,
              help = "Output file name prefix [default: sample name]"),

  make_option("--regions", type = "character", default = "all",
              help = "'all', or e.g. 'chr8,chr17:0-30000000' [default %default]"),
  make_option("--regions-file", type = "character", default = NULL,
              help = "BED file of regions (overrides --regions)"),
  make_option("--layout", type = "character", default = "separate",
              help = "separate | together | both  [default %default]"),

  make_option("--genes", type = "character", default = NULL,
              help = "Comma-separated HUGO gene symbols to annotate"),
  make_option("--title", type = "character", default = NULL,
              help = "Plot title [default: sample, purity/ploidy, region]"),
  make_option("--genome", type = "character", default = "hg38",
              help = "Genome build: hg38, hg19, T2T, mm10, mm39 [default %default]"),
  make_option("--max-cn", type = "double", default = 8,
              help = "Copy number axis ceiling [default %default]"),
  make_option("--format", type = "character", default = "pdf,png",
              help = "Comma-separated output formats [default %default]"),
  make_option("--width", type = "double", default = NULL, help = "Figure width (in)"),
  make_option("--height", type = "double", default = NULL, help = "Figure height (in)"),

  make_option("--baf-track", action = "store_true", default = FALSE,
              help = "Add a het-SNP BAF annotation panel (SAVANA only)"),
  make_option("--annotation-file", type = "character", default = NULL,
              help = "Explicit annotation source (het-SNP BED for --baf-track)"),
  make_option("--baf-max-points", type = "integer", default = 5000,
              help = "Het-SNPs drawn per panel after thinning [default %default]"),
  make_option("--extra", type = "character", default = NULL,
              help = paste("Escape hatch for any other ReConPlot() argument, e.g.",
                           "'size_text=6,curvature_intrachr_SVs=-0.2'")),
  make_option("--drop-na-minor", action = "store_true", default = FALSE,
              help = "Drop CN segments lacking a minor allele CN"),
  make_option("--keep-interchrom-strands", action = "store_true", default = FALSE,
              help = "Colour inter-chromosomal SVs by orientation instead of as TRA"),
  make_option("--write-tables", action = "store_true", default = FALSE,
              help = "Also write the harmonised CN/SV tables as TSV"),
  make_option("--seed", type = "integer", default = 1,
              help = "RNG seed (BAF thinning) [default %default]")
)

opt <- parse_args(OptionParser(
  usage = "%prog --source SOURCE --input DIR --outdir DIR [options]",
  option_list = option_list))

## optparse turns --cn-file into opt$`cn-file`; give everything snake_case names
## so parsers can use args$cn_file.
args <- opt
names(args) <- gsub("-", "_", names(args))

if (isTRUE(args$list_sources)) {
  print(list_parsers(), right = FALSE)
  quit(status = 0)
}
if (is.null(args$outdir)) stop("--outdir is required")
if (is.null(args$input) && is.null(args$cn_input) && is.null(args$cn_file)) {
  stop("Provide --input (a caller output directory) or explicit --cn-file/--sv-file")
}
set.seed(args$seed)

## --- 1. parse ---------------------------------------------------------------
parsed <- run_parsers(args)
meta <- parsed$meta %||% list()

## --- 2. harmonise / validate -----------------------------------------------
log_msg("validating inputs")
cn <- validate_cn(parsed$cn, drop_na_minor = isTRUE(args$drop_na_minor))
sv <- validate_sv(parsed$sv, interchrom_as_tra = !isTRUE(args$keep_interchrom_strands))
log_msg(sprintf("  %d CN segments, %d SV junctions retained", nrow(cn), nrow(sv)))
if (nrow(sv) > 0) {
  tab <- table(sv$strands)
  log_msg("  SV types: ", paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = "  "))
}

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
prefix <- args$prefix %||% meta$sample %||% args$source
if (isTRUE(args$write_tables)) {
  data.table::fwrite(cn, file.path(args$outdir, paste0(prefix, ".reconplot_cn.tsv")), sep = "\t")
  data.table::fwrite(sv, file.path(args$outdir, paste0(prefix, ".reconplot_sv.tsv")), sep = "\t")
  log_msg("  wrote harmonised tables to ", args$outdir)
}

## --- 3. regions -------------------------------------------------------------
regions <- if (!is.null(args$regions_file)) {
  read_regions_file(args$regions_file)
} else {
  parse_regions(args$regions, cn, genome = args$genome)
}
regions <- sanitize_regions(regions, cn)
log_msg(sprintf("%d region(s) to plot", nrow(regions)))

combined <- collapse_duplicate_chroms(regions)
region_sets <- switch(
  args$layout,
  separate = split(regions, seq_len(nrow(regions))),
  together = list(combined),
  both     = c(list(combined), split(regions, seq_len(nrow(regions)))),
  stop("--layout must be one of: separate, together, both")
)
region_sets <- lapply(region_sets, function(x) { rownames(x) <- NULL; x })

## --- 4. plot ----------------------------------------------------------------
genes <- if (!is.null(args$genes)) trimws(unlist(strsplit(args$genes, ","))) else NULL
formats <- trimws(unlist(strsplit(args$format, ",")))

## --extra lets any remaining ReConPlot() argument through without a dedicated flag
extra <- list(max.cn = args$max_cn, genome_version = args$genome)
if (!is.null(args$extra)) {
  for (kv in trimws(unlist(strsplit(args$extra, ",")))) {
    if (!nzchar(kv)) next
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("--extra entries must look like key=value: ", kv)
    val <- suppressWarnings(as.numeric(parts[2]))
    if (is.na(val)) {
      val <- switch(parts[2], "TRUE" = TRUE, "FALSE" = TRUE, parts[2])
      if (parts[2] == "FALSE") val <- FALSE
    }
    extra[[trimws(parts[1])]] <- val
  }
  log_msg("  extra ReConPlot args: ", paste(names(extra), unlist(extra), sep = "=", collapse = " "))
}

written <- plot_region_sets(
  region_sets, cn = cn, sv = sv, meta = meta,
  outdir = args$outdir, prefix = prefix,
  title = args$title, genes = genes,
  annotation_fn = parsed$annotation_fn,
  formats = formats, width = args$width, height = args$height,
  extra = extra
)

log_msg(sprintf("done: %d file(s) written to %s", length(written), args$outdir))
if (length(written) == 0) quit(status = 1)
