#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## run_shatterseek.R -- ShatterSeek chromothripsis calls from lrsomatic caller output.
##
## <caller output> --[ReConPlot wrapper parsers]--> {cn, sv} --[ShatterSeek format]--> shatterseek()
##   --> PCAWG high-confidence tiers (HC1/HC2/HC3) --> summary TSV + per-chromosome figures
##
## Callers are read by the parsers in assets/reconplot/R/parsers (--reconplot-dir), so ShatterSeek and
## ReConPlot see the same SVs and segments. ShatterSeek itself is upstream parklab/ShatterSeek, installed in the
## module's container with a build-time patch adding T2T-CHM13 (containers/shatterseek/).
## The tiers follow Cortes-Ciriano et al. 2020, Supplementary Note 1.2, as spelled out in parklab/ShatterSeek#22.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--reconplot-dir", type = "character", default = NULL,
              help = "ReConPlot wrapper directory whose R/ parsers read the caller output (required)"),
  make_option("--cn-source", type = "character", default = NULL, help = "Copy number parser: ascat, wakhan or savana"),
  make_option("--sv-source", type = "character", default = NULL, help = "SV parser: severus or savana"),
  make_option("--cn-input", type = "character", default = NULL, help = "Directory with the CN caller's files"),
  make_option("--sv-input", type = "character", default = NULL, help = "Directory with the SV caller's files [default: --cn-input]"),
  make_option("--sample", type = "character", default = NULL, help = "Sample name (also pins caller files)"),
  make_option("--prefix", type = "character", default = NULL, help = "Output file prefix [default: --sample]"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory [default %default]"),
  make_option("--genome", type = "character", default = "hg38", help = "hg38, chm13 or hg19 [default %default]"),
  make_option("--min-support", type = "double", default = 0, help = "Drop SVs with tumour read support below this [default %default]"),
  make_option("--min-svlen", type = "double", default = 0,
              help = "Drop intra-chromosomal SVs shorter than this many bp [default %default]"),
  make_option("--exclude-vntr", action = "store_true", default = FALSE, help = "Drop Severus SVs inside a VNTR"),
  make_option("--sv-format", type = "character", default = "auto", help = "SAVANA SV input: bedpe, vcf or auto [default %default]"),
  make_option("--ascat-cn-mode", type = "character", default = "segments",
              help = "ASCAT table: segments, raw or raw-fractional [default %default]"),
  make_option("--min-size", type = "integer", default = 1, help = "shatterseek() min.Size [default %default]"),
  make_option("--use", type = "character", default = "q", help = "Tier tests on BH q (q) or raw P values (p) [default %default]"),
  make_option("--alpha", type = "double", default = 0.2, help = "Significance level for the tiers [default %default]"),
  make_option("--seed", type = "integer", default = 1,
              help = "RNG seed (ShatterSeek's exponential tests draw random numbers) [default %default]"),
  make_option("--min-sv-hc1", type = "double", default = 6, help = "HC1 interleaved intra SVs [default %default]"),
  make_option("--min-osc-hc", type = "double", default = 7, help = "HC1/HC2 oscillating segments [default %default]"),
  make_option("--min-sv-hc2", type = "double", default = 3, help = "HC2 interleaved intra SVs [default %default]"),
  make_option("--min-tra-hc2", type = "double", default = 4, help = "HC2 translocations [default %default]"),
  make_option("--min-sv-hc3", type = "double", default = 40, help = "HC3 SVs incl. translocations [default %default]"),
  make_option("--no-hc3", action = "store_true", default = FALSE, help = "Skip the >= 40 SV tier (PCAWG: sarcomas only)"),
  make_option("--plot-all", action = "store_true", default = FALSE,
              help = "Draw every chromosome with a cluster, not only those with a tier")
)
args <- parse_args(OptionParser(usage = "%prog --reconplot-dir DIR --cn-source SRC --sv-source SRC --cn-input DIR [options]",
                                option_list = option_list))
names(args) <- gsub("-", "_", names(args))

if (is.null(args$reconplot_dir)) stop("--reconplot-dir is required")
for (f in c("utils_common.R", "utils_vcf.R", file.path("parsers", "registry.R"))) source(file.path(args$reconplot_dir, "R", f))
load_parsers(file.path(args$reconplot_dir, "R", "parsers"))
suppressPackageStartupMessages(library(ShatterSeek))

args$genome <- match.arg(args$genome, c("hg38", "chm13", "hg19"))
args$use <- match.arg(args$use, c("q", "p"))
args$sv_input <- args$sv_input %||% args$cn_input
prefix <- args$prefix %||% args$sample %||% "sample"
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

ANALYSED <- c(as.character(1:22), "X")             # ShatterSeek's chromosomes (Ensembl names)
STRAND_PAIRS <- c("+-", "-+", "++", "--")
SVTYPE_OF <- c("+-" = "DEL", "-+" = "DUP", "++" = "h2hINV", "--" = "t2tINV")

## --- 1. parse with the ReConPlot parsers -------------------------------------
parsed <- run_parsers(args)
cn <- validate_cn(parsed$cn)
sv <- validate_sv(parsed$sv, interchrom_as_tra = FALSE)   # keep each translocation's strands

## --- 2. SVs in ShatterSeek's format --------------------------------------------
to_shatterseek_sv <- function(sv) {
  if (nrow(sv) == 0) {
    return(data.frame(chrom1 = character(), pos1 = integer(), chrom2 = character(), pos2 = integer(),
                      SVtype = character(), strand1 = character(), strand2 = character()))
  }
  strands <- sv$strands
  by_class <- c(DEL = "+-", DUP = "-+", h2hINV = "++", t2tINV = "--")   # generic parser may give classes
  strands[strands %in% names(by_class)] <- by_class[strands[strands %in% names(by_class)]]
  keep <- strands %in% STRAND_PAIRS
  if (any(!keep)) {
    tab <- table(strands[!keep])
    log_msg("  dropped ", sum(!keep), " SVs without a two-breakend orientation (",
            paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "), ")")
  }
  out <- data.frame(chrom1 = sub("^chr", "", sv$chr1[keep]), pos1 = as.integer(sv$pos1[keep]),
                    chrom2 = sub("^chr", "", sv$chr2[keep]), pos2 = as.integer(sv$pos2[keep]),
                    strand1 = substr(strands[keep], 1, 1), strand2 = substr(strands[keep], 2, 2),
                    stringsAsFactors = FALSE)
  analysed <- out$chrom1 %in% ANALYSED & out$chrom2 %in% ANALYSED
  if (any(!analysed)) log_msg("  dropped ", sum(!analysed), " SVs with a breakend outside chromosomes 1-22 and X")
  out <- out[analysed, , drop = FALSE]
  # pos1 <= pos2, strands moving with their breakends (the order SVs() would impose)
  flip <- out$pos1 > out$pos2
  out[flip, c("chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2")] <-
    out[flip, c("chrom2", "pos2", "strand2", "chrom1", "pos1", "strand1")]
  out$SVtype <- ifelse(out$chrom1 == out$chrom2, unname(SVTYPE_OF[paste0(out$strand1, out$strand2)]), "TRA")
  out[order(match(out$chrom1, ANALYSED), out$pos1), c("chrom1", "pos1", "chrom2", "pos2", "SVtype", "strand1", "strand2")]
}

## --- 3. copy number in ShatterSeek's format -----------------------------------
to_shatterseek_cn <- function(cn) {
  d <- data.frame(chrom = sub("^chr", "", cn$chr), start = as.integer(cn$start), end = as.integer(cn$end),
                  total_cn = round(as.numeric(cn$copyNumber)), stringsAsFactors = FALSE)
  if (any(cn$copyNumber != d$total_cn, na.rm = TRUE)) {
    log_msg("  rounded ", sum(cn$copyNumber != d$total_cn, na.rm = TRUE), " non-integer total copy numbers")
  }
  d <- d[d$chrom %in% ANALYSED & !is.na(d$total_cn), , drop = FALSE]
  d <- d[order(match(d$chrom, ANALYSED), d$start, d$end), , drop = FALSE]
  out <- list()
  for (ch in unique(d$chrom)) {
    x <- d[d$chrom == ch, , drop = FALSE]
    # no overlaps: each segment starts after every earlier one ends (running maximum, so nested segments drop out)
    if (nrow(x) > 1) {
      prev_end <- c(-Inf, cummax(x$end)[-nrow(x)])
      x$start <- as.integer(pmax(x$start, prev_end + 1))
      x <- x[x$start <= x$end, , drop = FALSE]
    }
    # neighbours with equal copy number become one segment (ShatterSeek counts states row by row)
    run <- cumsum(c(TRUE, diff(x$total_cn) != 0))
    out[[ch]] <- data.frame(chrom = ch, start = tapply(x$start, run, min), end = tapply(x$end, run, max),
                            total_cn = tapply(x$total_cn, run, `[`, 1), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  if (is.null(out) || nrow(out) == 0) stop("No copy number segments on chromosomes 1-22 and X")
  out
}

sv_ss <- to_shatterseek_sv(sv)
cn_ss <- to_shatterseek_cn(cn)
log_msg(sprintf("ShatterSeek input: %d SVs (%s), %d CN segments after merging equal neighbours",
                nrow(sv_ss), paste(sprintf("%s=%d", names(table(sv_ss$SVtype)), as.integer(table(sv_ss$SVtype))), collapse = " "),
                nrow(cn_ss)))
fwrite(sv_ss, file.path(args$outdir, paste0(prefix, ".shatterseek_sv_input.tsv")), sep = "\t")
fwrite(cn_ss, file.path(args$outdir, paste0(prefix, ".shatterseek_cn_input.tsv")), sep = "\t")

## --- 4. ShatterSeek + PCAWG tiers ---------------------------------------------
## PCAWG high-confidence tiers (Cortes-Ciriano et al. 2020, Supplementary Note 1.2; parklab/ShatterSeek#22). Tests
## are significant at <= alpha on BH q (across the sample's chromosomes) or raw P; NA counts as a fail.
classify_pcawg <- function(d, use, alpha, min_sv_hc1, min_osc_hc, min_sv_hc2, min_tra_hc2, min_sv_hc3, hc3) {
  for (stem in c("fragment_joins", "chr_breakpoint_enrichment", "exp_chr", "exp_cluster")) {
    p_col <- if (stem == "chr_breakpoint_enrichment") stem else paste0("pval_", stem)
    d[[paste0("q_", stem)]] <- p.adjust(d[[p_col]], method = "BH")
  }
  test <- function(stem) d[[if (use == "q") paste0("q_", stem) else if (stem == "chr_breakpoint_enrichment") stem
                            else paste0("pval_", stem)]]
  ok <- function(x) !is.na(x) & x

  n_intra <- d$number_DEL + d$number_DUP + d$number_h2hINV + d$number_t2tINV
  osc <- d$max_number_oscillating_CN_segments_2_states
  joins_random <- ok(test("fragment_joins") > alpha)
  clustered <- ok(test("chr_breakpoint_enrichment") <= alpha) | ok(test("exp_chr") <= alpha)

  d$HC1 <- ok(n_intra >= min_sv_hc1) & ok(osc >= min_osc_hc) & joins_random & clustered
  d$HC2 <- ok(n_intra >= min_sv_hc2) & ok(d$number_TRA >= min_tra_hc2) & ok(osc >= min_osc_hc) & joins_random
  d$HC3 <- hc3 & ok(d$clusterSize_including_TRA >= min_sv_hc3) & joins_random

  tiers <- c("HC1", "HC2", "HC3")
  d$chromothripsis_criteria <- apply(as.matrix(d[, tiers]), 1, function(r) paste(tiers[r], collapse = ","))
  d$chromothripsis_confidence <- ifelse(d$HC1 | d$HC2 | d$HC3, "high", "none")
  d
}

set.seed(args$seed)
res <- shatterseek(
  SV.sample = SVs(chrom1 = sv_ss$chrom1, pos1 = as.numeric(sv_ss$pos1), chrom2 = sv_ss$chrom2, pos2 = as.numeric(sv_ss$pos2),
                  SVtype = sv_ss$SVtype, strand1 = sv_ss$strand1, strand2 = sv_ss$strand2),
  seg.sample = CNVsegs(chrom = cn_ss$chrom, start = as.numeric(cn_ss$start), end = as.numeric(cn_ss$end),
                       total_cn = as.numeric(cn_ss$total_cn)),
  min.Size = args$min_size, genome = args$genome)
res@chromSummary <- classify_pcawg(
  res@chromSummary, use = args$use, alpha = args$alpha, min_sv_hc1 = args$min_sv_hc1, min_osc_hc = args$min_osc_hc,
  min_sv_hc2 = args$min_sv_hc2, min_tra_hc2 = args$min_tra_hc2, min_sv_hc3 = args$min_sv_hc3, hc3 = !isTRUE(args$no_hc3))

summary <- res@chromSummary
summary$inter_other_chroms_coords_all <- trimws(gsub("[\n;]+", ";", summary$inter_other_chroms_coords_all))
out <- cbind(data.frame(sample = args$sample %||% prefix, cn_source = args$cn_source, sv_source = args$sv_source,
                        stringsAsFactors = FALSE), summary)
out$chrom <- paste0("chr", out$chrom)
out$start <- as.integer(out$start)
out$end <- as.integer(out$end)
fwrite(out, file.path(args$outdir, paste0(prefix, ".shatterseek.tsv")), sep = "\t", na = "NA", quote = FALSE)
called <- summary[summary$chromothripsis_confidence != "none", , drop = FALSE]
log_msg(sprintf("tiers: %s", if (nrow(called)) paste(sprintf("chr%s %s (%s)", called$chrom, called$chromothripsis_confidence,
                                                                called$chromothripsis_criteria), collapse = "; ") else "none"))

## --- 5. figures -----------------------------------------------------------------
draw <- summary$chrom[if (isTRUE(args$plot_all)) !is.na(summary$start) else summary$chromothripsis_confidence != "none"]
if (length(draw)) dir.create(file.path(args$outdir, "plots"), showWarnings = FALSE)
for (ch in draw) {
  f <- file.path(args$outdir, "plots", sprintf("%s_chr%s.pdf", prefix, ch))
  ok <- tryCatch({
    p <- plot_chromothripsis(ShatterSeek_output = res, chr = ch, sample_name = prefix, genome = args$genome)
    pdf(f, width = 9, height = 11)
    gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 4, ncol = 1, heights = c(0.2, 0.4, 0.4, 0.4))
    dev.off()
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    unlink(f)
    log_msg("ERROR: figure for chr", ch, " not drawn: ", conditionMessage(e))
    FALSE
  })
  if (ok) log_msg("  drew ", basename(f))
}
log_msg("done")
