#!/usr/bin/env Rscript
# Build-time patch adding T2T-CHM13 to upstream parklab/ShatterSeek. Run from the unpacked package root:
#   Rscript patch_chm13.R <chm13v2.0_cytobands_allchrs.bed>
#
# Upstream only knows hg19 and hg38, through three genome checks and two per-genome tables in R/sysdata.rda:
#   info_mappa_hg38 (chromosome sizes for the breakpoint-enrichment test) and chr_info_hg38 (ideogram bands).
# This adds info_mappa_chm13 (CHM13 has no gaps, so the sizes are the chromosome lengths) and chr_info_chm13 in the same
# layout, and lets the three checks accept "chm13". The statistics are untouched. Fails if any edit does not apply.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, file.exists(args[1]))

cyto <- read.delim(args[1], header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("seqnames", "start", "end", "name", "gieStain"))

env <- new.env()
load("R/sysdata.rda", envir = env)
stopifnot(exists("info_mappa_hg38", envir = env), exists("chr_info_hg38", envir = env))

# chr_info_chm13: the bands, in chr_info_hg38's columns and chromosome naming
hg38 <- env$chr_info_hg38
cyto <- cyto[sub("^chr", "", cyto$seqnames) %in% sub("^chr", "", unique(as.character(hg38$seqnames))), ]
if (!any(grepl("^chr", hg38$seqnames))) cyto$seqnames <- sub("^chr", "", cyto$seqnames)
cyto$width <- cyto$end - cyto$start + 1L   # as in chr_info_hg38
cyto$strand <- "*"
cyto$y <- 1L
chr_info_chm13 <- cyto[, intersect(names(hg38), names(cyto)), drop = FALSE]
stopifnot(identical(names(chr_info_chm13), names(hg38)))
if (is.factor(hg38$seqnames)) chr_info_chm13$seqnames <- factor(chr_info_chm13$seqnames, levels = unique(chr_info_chm13$seqnames))

# info_mappa_chm13: one row per chromosome of info_mappa_hg38 (1-22, X, Y, M), tot = chromosome length
mappa_hg38 <- env$info_mappa_hg38
chm13_len <- c(tapply(cyto$end, sub("^chr", "", cyto$seqnames), max), M = 16569)   # chrM (rCRS) has no bands
info_mappa_chm13 <- mappa_hg38
info_mappa_chm13$tot <- as.numeric(chm13_len[sub("^chr", "", as.character(mappa_hg38$V1))])
stopifnot(!anyNA(info_mappa_chm13$tot), all(info_mappa_chm13$tot > 0))
for (col in setdiff(names(info_mappa_chm13), c("V1", "tot"))) info_mappa_chm13[[col]] <- NA

assign("chr_info_chm13", chr_info_chm13, envir = env)
assign("info_mappa_chm13", info_mappa_chm13, envir = env)
save(list = ls(env, all.names = TRUE), envir = env, file = "R/sysdata.rda", compress = "bzip2")

# the three genome checks: (file, exact upstream text, replacement)
edits <- list(
  list("R/shatterseek.R",
       'if ( !(as.character(genome) %in% c("hg19","hg38"))){stop("Reference genome assembly is not supported (Use hg19 or hg38)")}',
       'if ( !(as.character(genome) %in% c("hg19","hg38","chm13"))){stop("Reference genome assembly is not supported (Use hg19, hg38 or chm13)")}'),
  list("R/function_criteria.R",
       'if(genome == "hg38"){',
       'if(genome == "chm13"){\n              prob_cand = info_mappa_chm13$tot[which(info_mappa_chm13$V1 == cand)] / sum(as.numeric(info_mappa_chm13$tot))\n            }else if(genome == "hg38"){'),
  list("R/plot_function.R",
       '}else if(genome == "hg38"){\n\t  chr_info <- chr_info_hg38',
       '}else if(genome == "chm13"){\n\t  chr_info <- chr_info_chm13\n\t}else if(genome == "hg38"){\n\t  chr_info <- chr_info_hg38')
)
for (e in edits) {
  src <- paste(readLines(e[[1]], warn = FALSE), collapse = "\n")
  n <- lengths(regmatches(src, gregexpr(e[[2]], src, fixed = TRUE)))
  if (n != 1) stop(sprintf("%s: expected the upstream text once, found it %d times", e[[1]], n))
  writeLines(sub(e[[2]], e[[3]], src, fixed = TRUE), e[[1]])
}
cat("ShatterSeek patched for chm13:", nrow(chr_info_chm13), "bands,", nrow(info_mappa_chm13), "chromosome sizes\n")
