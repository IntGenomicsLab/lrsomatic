suppressPackageStartupMessages({
  library(data.table)
})

# Parse Severus somatic VCF for circos plot data.
# Returns list: $translocations (BND records) and $nontrans (DEL/DUP/INV/INS)
parse_severus_vcf = function(vcf_file) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
    return(list(translocations = data.table(), nontrans = data.table()))
  }

  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  # fread() errors (rather than returning 0 rows) when skip lands exactly on
  # the last line of the file, i.e. a VCF with no variant records at all.
  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
          select = 1:8,
          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0) {
    return(list(translocations = data.table(), nontrans = data.table()))
  }

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract INFO sub-fields
  .info_val = function(info_vec, key) {
    pattern = paste0("(?:^|;)", key, "=([^;]+)")
    m = regmatches(info_vec, regexpr(pattern, info_vec, perl = TRUE))
    ifelse(nchar(m) > 0, sub(paste0(".*="), "", m), NA_character_)
  }

  dt[, SVTYPE := .info_val(INFO, "SVTYPE")]
  dt[grepl("END=",   INFO, fixed = TRUE), END   := as.integer(.info_val(INFO[grepl("END=",   INFO, fixed = TRUE)], "END"))]
  dt[grepl("SVLEN=", INFO, fixed = TRUE), SVLEN := as.integer(.info_val(INFO[grepl("SVLEN=", INFO, fixed = TRUE)], "SVLEN"))]

  # BND partner chromosome/position from ALT field
  # ALT format examples: "N[chr7:24547089[" or "]chr7:24547089]N"
  dt[SVTYPE == "BND", CHROM2 := {
    m = regmatches(ALT, regexpr("chr[^:]+", ALT, perl = TRUE))
    ifelse(nchar(m) > 0, m, NA_character_)
  }]
  dt[SVTYPE == "BND", POS2 := as.integer(regmatches(ALT, regexpr("(?<=:)\\d+", ALT, perl = TRUE)))]

  # Insertions have no END — use POS
  dt[SVTYPE == "INS" | is.na(END), END := POS]

  # Colours and y-positions for non-BND SV track
  SV_COL  = c(INS = "#f97e02", DEL = "#020272", INV = "#e7cc02", DUP = "#e41a1c")
  SV_YPOS = c(INS = 1.0,       DEL = 0.66,      INV = 0.33,      DUP = 0.05)
  dt[SVTYPE %in% names(SV_COL),  circos_col := SV_COL[SVTYPE]]
  dt[SVTYPE %in% names(SV_YPOS), circos_pos := SV_YPOS[SVTYPE]]

  translocations = dt[SVTYPE == "BND" & !is.na(CHROM2) & !is.na(POS2),
    .(chrom = CHROM, pos = POS, chrom2 = CHROM2, pos2 = POS2)]

  nontrans = dt[SVTYPE != "BND",
    .(chrom = CHROM, pos = POS, end = END, svtype = SVTYPE,
      svlen = SVLEN, circos_pos, circos_col)]

  list(translocations = translocations, nontrans = nontrans)
}

# Parse the somatic Severus VCF into one row per SV (id, svtype, coords, length, VAF).
# Used as the input to build_sv_table_from_vep() — a lighter-weight companion to
# parse_severus_vcf() above, which instead returns circos-ready translocation/non-BND tracks.
parse_severus_somatic_records = function(vcf_file) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(data.table())

  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE, select = 1:10,
          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                        "FORMAT", "SAMPLE1")),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0) return(data.table())

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  .info_val = function(info_vec, key) {
    pattern = paste0("(?:^|;)", key, "=([^;]+)")
    m = regmatches(info_vec, regexpr(pattern, info_vec, perl = TRUE))
    ifelse(nchar(m) > 0, sub(paste0(".*="), "", m), NA_character_)
  }

  dt[, SVTYPE := .info_val(INFO, "SVTYPE")]
  dt[grepl("END=",   INFO, fixed = TRUE), END   := as.integer(.info_val(INFO[grepl("END=",   INFO, fixed = TRUE)], "END"))]
  dt[grepl("SVLEN=", INFO, fixed = TRUE), SVLEN := as.integer(.info_val(INFO[grepl("SVLEN=", INFO, fixed = TRUE)], "SVLEN"))]

  # Insertions and BNDs have no END — use POS
  dt[SVTYPE %in% c("INS", "BND") | is.na(END), END := POS]

  # VAF from FORMAT/SAMPLE1 (format string is uniform for Severus output, but split
  # format-group by format-group defensively, as in parse_caller_vcf())
  fmt_groups = unique(dt$FORMAT)
  vaf_list = vector("numeric", nrow(dt))
  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields = strsplit(fmt, ":", fixed = TRUE)[[1]]
    vaf_idx = match("VAF", fields)
    split_s = strsplit(dt$SAMPLE1[idx_rows], ":", fixed = TRUE)
    vaf_list[idx_rows] = if (!is.na(vaf_idx)) {
      vapply(split_s, function(x)
        if (length(x) >= vaf_idx) suppressWarnings(as.numeric(x[vaf_idx])) else NA_real_,
        numeric(1))
    } else NA_real_
  }
  dt[, VAF := vaf_list]

  dt[, .(id = ID, svtype = SVTYPE, chrom = CHROM, start = POS, end = END,
         sv_len = SVLEN, vaf = VAF)]
}

# Build the SV display table by joining VEP CSQ gene annotations (from the SV VEP VCF,
# `parse_vep_vcf()` from R/parse_smallvariants.R) onto the somatic Severus SVs by locus.
# This is the primary path when a VEP SV VCF is available; build_sv_table() below (fed by
# the gene-annotated TSV) is the fallback for pipelines that don't produce a VEP SV VCF.
build_sv_table_from_vep = function(somatic_vcf, vep_sv_vcf) {
  somatic = parse_severus_somatic_records(somatic_vcf)
  if (nrow(somatic) == 0) return(data.table())

  vep = parse_vep_vcf(vep_sv_vcf)
  if (is.null(vep) || nrow(vep) == 0) {
    somatic[, `:=`(gene_hits = NA_character_, consequence = NA_character_, impact = NA_character_)]
    return(somatic[, .(id, gene_hits, svtype, chrom, start, end, sv_len, vaf, consequence, impact)])
  }

  # Keep the highest-impact annotation per locus, and collapse all distinct gene symbols
  # hit at that locus into one comma-joined column.
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep[, impact_rank := impact_rank[impact]]
  vep[is.na(impact_rank), impact_rank := 5L]
  setorder(vep, impact_rank)

  agg = vep[, .(
    gene_hits   = paste(unique(symbol[!is.na(symbol) & nzchar(symbol)]), collapse = ","),
    consequence = consequence[1],
    impact      = impact[1]
  ), by = .(chrom, pos)]
  agg[!nzchar(gene_hits), gene_hits := NA_character_]

  merged = merge(somatic, agg, by.x = c("chrom", "start"), by.y = c("chrom", "pos"), all.x = TRUE)
  merged[, .(id, gene_hits, svtype, chrom, start, end, sv_len, vaf, consequence, impact)]
}

# Parse the gene-annotated Severus TSV (filtered_SV2/SV_filtered_with_gene_annotations.tsv)
parse_severus_gene_tsv = function(tsv_file) {
  if (is.null(tsv_file) || !file.exists(tsv_file)) return(NULL)
  dt = fread(tsv_file, sep = "\t", header = TRUE, fill = TRUE)
  setnames(dt, toupper(names(dt)))

  if ("START_CHROM" %in% names(dt)) dt[, START_CHROM := ensure_chr_prefix(START_CHROM)]
  if ("END_CHROM"   %in% names(dt)) dt[, END_CHROM   := ensure_chr_prefix(END_CHROM)]

  # Gene column: prefer NHL hits
  gene_col = if ("NHL_GENE_HITS"    %in% names(dt)) "NHL_GENE_HITS"
             else if ("COSMIC_GENE_HITS" %in% names(dt)) "COSMIC_GENE_HITS"
             else NULL
  dt[, gene_hits := if (!is.null(gene_col)) get(gene_col) else NA_character_]
  dt
}

# Build the SV display table.
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all SVs (one row each).
build_sv_table = function(sv_tsv, gene_panel = NULL) {
  if (is.null(sv_tsv) || nrow(sv_tsv) == 0) return(data.table())

  sv_tsv = copy(sv_tsv)

  # When no panel is supplied, return one row per SV without explosion
  if (is.null(gene_panel)) {
    display_cols = intersect(
      c("ID", "SVTYPE", "DETAILED_TYPE",
        "START_CHROM", "START_POS", "END_CHROM", "END_POS",
        "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
        "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
      names(sv_tsv)
    )
    return(sv_tsv[, ..display_cols])
  }

  # Panel-filtered path: explode multi-gene gene_hits, filter, return one row per gene×SV
  sv_tsv[, .ridx := .I]

  sv_long = sv_tsv[, {
    raw = as.character(gene_hits[1])
    genes = unique(trimws(unlist(strsplit(raw, "[;,]+"))))
    genes = genes[nchar(genes) > 0 & genes != "-" & toupper(genes) != "NA"]
    if (length(genes) == 0) genes = NA_character_
    list(gene = genes)
  }, by = .ridx]

  sv_long = merge(sv_long, sv_tsv, by = ".ridx")
  sv_long[, .ridx := NULL]
  sv_tsv[,  .ridx := NULL]

  # Filter by gene panel (always, even if panel is empty)
  sv_long = sv_long[!is.na(gene) & gene %in% gene_panel]
  if (nrow(sv_long) == 0) return(data.table())

  display_cols = intersect(
    c("gene", "ID", "SVTYPE", "DETAILED_TYPE",
      "START_CHROM", "START_POS", "END_CHROM", "END_POS",
      "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
      "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
    names(sv_long)
  )
  sv_long[, ..display_cols]
}
