suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Values of LRSomatic's INFO/CALLER tag that denote a somatic caller. The pipeline's
# *_SOMATIC_VEP.vcf.gz is a merged multi-caller VCF in which germline callers
# (deepvariant, clair3) supply the overwhelming majority of records, so the somatic
# table has to be filtered on this tag. ClairS is tagged "clairs" in matched mode and
# "clairs-to" in tumour-only mode.
SOMATIC_CALLERS = c("clairs", "clairs-to", "clairsto", "deepsomatic")

# Derive dbsnp/cosmic columns from a VEP "Existing_variation" column (semicolon- or
# comma-joined list of IDs, e.g. "rs123&COSV456"). Shared by parse_vep_text/parse_vep_vcf.
derive_dbsnp_cosmic = function(dt) {
  dt[, dbsnp := sub("(rs[0-9]+).*", "\\1", existing)]
  dt[!grepl("^rs", dbsnp, perl = TRUE), dbsnp := NA_character_]

  dt[, cosmic := sub(".*(COS[VM][0-9]+).*", "\\1", existing)]
  dt[!grepl("^COS", cosmic, perl = TRUE), cosmic := NA_character_]
  dt
}

# Dispatch to the right VEP parser based on actual file contents — both forms ship as
# "*_SOMATIC_VEP.vcf.gz" so the filename alone doesn't tell you which one you have.
# - VEP default text output: "##"-commented header, column line starts with "#Uploaded_variation"
# - genuine VCF w/ CSQ INFO field: "##fileformat=VCFv4.2", column line starts with "#CHROM"
parse_vep = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  con = gzfile(vep_file, "rb")
  is_vcf = FALSE
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    if (startsWith(line, "#CHROM")) { is_vcf = TRUE; break }
  }
  close(con)

  if (is_vcf) parse_vep_vcf(vep_file) else parse_vep_text(vep_file)
}

# Parse the VEP default text output (tab-delimited, ##-commented header, NOT a VCF).
# Returns a data.table with one row per consequence per variant.
parse_vep_text = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Count meta-lines (start with ##) to find the column-header line
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    skip_n = skip_n + 1L
  }
  close(con)

  dt = tryCatch(
    fread(vep_file, skip = skip_n, sep = "\t", header = TRUE,
          col.names = function(x) gsub("^#", "", x)),
    error = function(e) {
      message("Failed to parse VEP file: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  setnames(dt, old = "Uploaded_variation", new = "variant_id", skip_absent = TRUE)
  setnames(dt, old = "Gene",              new = "gene_id",     skip_absent = TRUE)
  setnames(dt, old = "Consequence",       new = "consequence", skip_absent = TRUE)

  # Coordinates and alleles both come from variant_id ("chr1_3506_A/G") wherever it has
  # that canonical shape, which is what VEP synthesises for VCF input without an ID.
  # Mixing the two sources is not safe: for an insertion reported in VEP's dash form,
  # Location's start is one base left of the position variant_id names
  # ("chr1_197488_-/G" has Location "chr1:197487-197488"), which then fails to join to
  # anything. Location remains the fallback for rows carrying a real VCF ID instead.
  vid = "^.+_[0-9]+_[^_]+/[^_]+$"
  dt[, from_vid := grepl(vid, variant_id)]

  dt[, chrom := ifelse(from_vid, sub("_[0-9]+_[^_]+$", "", variant_id),
                                 sub(":.*", "", Location))]
  dt[, pos   := as.integer(ifelse(from_vid, sub(".*_([0-9]+)_[^_]+$", "\\1", variant_id),
                                            sub(".*:(\\d+).*", "\\1", Location)))]
  dt[, chrom := ensure_chr_prefix(chrom)]

  dt[, ref := sub(".*_([^/]+)/.*", "\\1", variant_id)]
  dt[, alt := sub(".*/", "", variant_id)]

  # Parse VEP Extra key=value field
  dt[, symbol   := extract_extra_key(Extra, "SYMBOL")]
  dt[, impact   := extract_extra_key(Extra, "IMPACT")]
  dt[, existing := extract_extra_key(Extra, "Existing_variation")]
  dt[, sift     := extract_extra_key(Extra, "SIFT")]
  dt[, polyphen := extract_extra_key(Extra, "PolyPhen")]
  dt[, hgvsp    := extract_extra_key(Extra, "HGVSp")]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt = derive_dbsnp_cosmic(dt)

  # No per-variant caller in the text format; keep the column for contract parity
  # with parse_vep_vcf().
  dt[, caller := NA_character_]

  # Record identity, for contract parity with parse_vep_vcf()'s `id`. The text format
  # has no VCF ID column, so this is VEP's own Uploaded_variation name — never a
  # caller-assigned ID, which is why the SV table's ID-keyed join only ever sees the
  # VCF path.
  dt[, id := variant_id]

  # chrom/pos/ref/alt here come from `variant_id`, i.e. VEP's own notation — indels
  # already shifted one base right and possibly dash-form. See coord_space in
  # build_variant_table().
  dt[, coord_space := "vep"]

  dt
}

# Parse a genuine VCF carrying VEP annotation in a CSQ INFO field (VEP run with --vcf,
# as opposed to the default text output handled by parse_vep_text()).
# Returns a data.table with one row per gene/transcript annotation per variant, using the
# same column contract as parse_vep_text(): chrom, pos, ref, alt, symbol, gene_id,
# consequence, impact, hgvsp, existing, dbsnp, cosmic, sift, polyphen, caller.
parse_vep_vcf = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Skip header to #CHROM, capturing the CSQ field order from its INFO meta-line
  # (e.g. "...Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...") and noting whether
  # the file carries a per-record CALLER tag.
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  csq_format = NULL
  has_caller_info = FALSE
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##INFO=<ID=CSQ")) {
      m = regmatches(line, regexpr("Format: [^\"]+", line))
      if (length(m) > 0) csq_format = strsplit(sub("^Format: ", "", m), "|", fixed = TRUE)[[1]]
    }
    if (startsWith(line, "##INFO=<ID=CALLER")) has_caller_info = TRUE
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  if (is.null(csq_format)) {
    message("Failed to parse VEP VCF: no CSQ Format found in header")
    return(NULL)
  }

  dt = tryCatch(
    fread(vep_file, skip = skip_n + 1L, sep = "\t", header = FALSE, select = 1:8,
          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")),
    error = function(e) {
      message("Failed to parse VEP VCF: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Drop records the callers themselves rejected (RefCall / LowQual / GERMLINE).
  dt = dt[FILTER %in% c("PASS", ".")]
  if (nrow(dt) == 0) return(NULL)

  # Fixed-string splits (no regex backtracking) — much faster than sub() on the full
  # semicolon-delimited INFO field across millions of rows.
  dt[, CSQ := tstrsplit(INFO, "CSQ=", fixed = TRUE, keep = 2L)[[1]]]
  dt[, CSQ := tstrsplit(CSQ,  ";",    fixed = TRUE, keep = 1L)[[1]]]

  # Keep only somatic-caller records. Older single-caller VEP VCFs carry no CALLER tag
  # at all, and filtering those on it would empty the table, so it is skipped for them.
  if (has_caller_info) {
    dt[, caller := tstrsplit(INFO, "CALLER=", fixed = TRUE, keep = 2L)[[1]]]
    dt[, caller := tstrsplit(caller, ";", fixed = TRUE, keep = 1L)[[1]]]
    dt = dt[caller %in% SOMATIC_CALLERS]
    if (nrow(dt) == 0) return(NULL)
  } else {
    dt[, caller := NA_character_]
  }

  # A variant reported by two somatic callers appears as two records; collapse to one
  # row per variant so the CSQ expansion below doesn't duplicate every annotation.
  # ID rides along: it is what links an SV annotation record back to the caller record
  # it came from (see build_sv_table_from_vep()), and is more reliable than the locus.
  dt = dt[, .(CSQ = CSQ[1], id = ID[1],
              caller = paste(sort(unique(caller)), collapse = ",")),
          by = .(CHROM, POS, REF, ALT)]
  dt[id == ".", id := NA_character_]

  # One row per gene/transcript annotation (comma-separated CSQ entries)
  dt_long = dt[, .(csq_entry = unlist(strsplit(CSQ, ",", fixed = TRUE))),
               by = .(CHROM, POS, REF, ALT, id, caller)]
  if (nrow(dt_long) == 0) return(NULL)

  # Split each entry on "|", materialising only the fields actually used below
  # (tstrsplit returns columns in the order requested via `keep`).
  need = c("Consequence", "IMPACT", "SYMBOL", "Gene", "HGVSp",
           "Existing_variation", "SIFT", "PolyPhen")
  keep_idx = match(need, csq_format)
  ok = !is.na(keep_idx)
  parts = tstrsplit(dt_long$csq_entry, "|", fixed = TRUE, fill = NA_character_,
                    keep = keep_idx[ok])
  names(parts) = need[ok]

  get_field = function(name) {
    if (is.null(parts[[name]])) rep(NA_character_, nrow(dt_long)) else parts[[name]]
  }

  dt_long[, consequence := gsub("&", ",", get_field("Consequence"))]
  dt_long[, impact      := get_field("IMPACT")]
  dt_long[, symbol      := get_field("SYMBOL")]
  dt_long[, gene_id     := get_field("Gene")]
  dt_long[, hgvsp       := get_field("HGVSp")]
  dt_long[, existing    := get_field("Existing_variation")]
  dt_long[, sift        := get_field("SIFT")]
  dt_long[, polyphen    := get_field("PolyPhen")]

  # Blank fields are "" in CSQ, not NA — normalise for consistency with parse_vep_text().
  # `caller` collapses to "" when the file had no CALLER tag, so it normalises here too.
  for (col in c("symbol", "impact", "hgvsp", "existing", "gene_id", "sift", "polyphen",
                "caller")) {
    dt_long[get(col) == "", (col) := NA_character_]
  }

  dt_long[, chrom := ensure_chr_prefix(CHROM)]
  dt_long[, pos   := as.integer(POS)]
  dt_long[, ref   := REF]
  dt_long[, alt   := ALT]

  # Unlike parse_vep_text(), these are the VCF's own POS/REF/ALT — the *original* input
  # coordinate, not VEP's shifted one. Keying them as VEP-space is what stopped every
  # indel on this path from joining. See build_variant_table().
  dt_long[, coord_space := "vcf"]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt_long = derive_dbsnp_cosmic(dt_long)

  dt_long[, .(chrom, pos, ref, alt, id, symbol, gene_id, consequence, impact, hgvsp,
              existing, dbsnp, cosmic, sift, polyphen, caller, coord_space)]
}

# Which sample column of a VCF to read, given its sample names.
# One sample is unambiguous; otherwise prefer the one named for this report, then the
# first that is not obviously the normal. Falling through to the first column is a guess,
# so say so — silently reporting the normal's VAF is the failure this exists to prevent.
pick_sample_column = function(samples, sample_id = NULL, vcf_file = "") {
  if (length(samples) == 1L) return(1L)

  if (!is.null(sample_id) && nzchar(sample_id)) {
    i = match(sample_id, samples)
    if (!is.na(i)) return(i)
  }
  not_normal = which(!grepl("^(normal|blood|germline|ref)", samples, ignore.case = TRUE))
  if (length(not_normal) > 0) {
    if (length(not_normal) > 1L || is.null(sample_id))
      warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
              "]; reading '", samples[not_normal[1]], "'.")
    return(not_normal[1])
  }
  warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
          "] and none looks like a tumour; reading '", samples[1], "'.")
  1L
}

# The allele fraction, from whichever FORMAT tag this VCF happens to use.
#
# LRSomatic *renames* this tag according to which caller was prioritised — see
# STANDARDIZE_AF in subworkflows/local/small_variant_consensus.nf: AF -> VAF for
# deepvariant/deepsomatic, VAF -> AF for clair — and a `consensus` merge skips the
# renaming entirely, so both names reach this parser. Hard-coding "AF" is what produced
# an entirely empty VAF column on a DeepSomatic run. AD is the last resort: every caller
# in the pipeline emits it, and it is the only recoverable source when neither tag exists.
allele_fraction = function(field, fields) {
  for (tag in c("AF", "VAF")) {
    if (tag %in% fields) {
      # Multi-allelic records carry one value per ALT; the first is the one that pairs
      # with the record's first ALT allele, which is all this table joins on.
      v = suppressWarnings(as.numeric(sub(",.*$", "", field(tag))))
      if (any(!is.na(v))) return(v)
    }
  }
  if ("AD" %in% fields) {
    ad = strsplit(field("AD"), ",", fixed = TRUE)
    return(vapply(ad, function(x) {
      n = suppressWarnings(as.numeric(x))
      if (length(n) < 2L || anyNA(n[1:2]) || sum(n[1:2]) == 0) return(NA_real_)
      n[2] / sum(n[1:2])
    }, numeric(1)))
  }
  rep(NA_real_, length(field("GT")))
}

# Parse a raw caller VCF for variant coordinates, VAF, depth and phasing.
# `vcf_file` may be several paths, in which case they are read and stacked into one
# table — ClairS splits its matched-mode output into snvs.vcf.gz + indel.vcf.gz.
# Returns data.table with: chrom, pos, ref, alt, vaf, dp, gt, ps, caller
#
# `sample_id` picks the sample column by name. Reading it positionally reports the
# *normal's* VAF/DP on a matched two-sample VCF — wrong numbers, no error, and nothing in
# the report that would reveal it.
parse_caller_vcf = function(vcf_file, caller_name = "unknown", sample_id = NULL) {
  if (is.null(vcf_file)) return(NULL)
  if (length(vcf_file) > 1) {
    parts = lapply(vcf_file, parse_caller_vcf, caller_name = caller_name,
                   sample_id = sample_id)
    parts = parts[!vapply(parts, is.null, logical(1))]
    return(if (length(parts) > 0) rbindlist(parts) else NULL)
  }
  if (!file.exists(vcf_file)) return(NULL)

  # Count header lines, keeping the #CHROM line itself — it carries the sample names.
  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  chrom_line = NULL
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) { chrom_line = line; break }
    skip_n = skip_n + 1L
  }
  close(con)
  if (is.null(chrom_line)) return(NULL)

  header  = strsplit(chrom_line, "\t", fixed = TRUE)[[1]]
  samples = if (length(header) > 9L) header[10:length(header)] else character(0)
  # A sites-only VCF has no genotypes to read. Return NULL rather than letting fread()
  # error on `select = 1:10` — every other parser here degrades gracefully.
  if (length(samples) == 0) return(NULL)

  sample_col = pick_sample_column(samples, sample_id, vcf_file)

  # 1:9 are the fixed VCF columns; the chosen sample sits at 9 + its index.
  col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                "SAMPLE1")
  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
          select = c(1:9, 9L + sample_col), col.names = col_names),
    error = function(e) {
      warning("Could not read VCF '", vcf_file, "': ", conditionMessage(e))
      NULL
    })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract the allele fraction, DP, GT and PS from the FORMAT + SAMPLE1 columns.
  # Work format-group by format-group to avoid splitting every single row redundantly.
  fmt_groups = unique(dt$FORMAT)
  vaf_list = rep(NA_real_,      nrow(dt))
  dp_list  = rep(NA_integer_,   nrow(dt))
  gt_list  = rep(NA_character_, nrow(dt))
  ps_list  = rep(NA_character_, nrow(dt))

  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields  = strsplit(fmt, ":", fixed = TRUE)[[1]]
    split_s = strsplit(dt$SAMPLE1[idx_rows], ":", fixed = TRUE)

    # One FORMAT field, by name, across this group's rows
    field = function(name) {
      i = match(name, fields)
      if (is.na(i)) return(rep(NA_character_, length(split_s)))
      vapply(split_s, function(x) if (length(x) >= i) x[i] else NA_character_,
             character(1))
    }

    vaf_list[idx_rows] = allele_fraction(field, fields)
    dp_list[idx_rows]  = suppressWarnings(as.integer(field("DP")))
    gt_list[idx_rows]  = field("GT")
    ps_list[idx_rows]  = field("PS")
  }

  # Unphased records carry "." for PS and a "/"-separated GT. Blank the placeholders so
  # the report shows an empty cell rather than a bare ".".
  ps_list[!is.na(ps_list) & ps_list == "."] = NA_character_
  gt_list[!is.na(gt_list) & gt_list %in% c(".", "./.")] = NA_character_

  data.table(chrom = dt$CHROM, pos = dt$POS, ref = dt$REF, alt = dt$ALT,
             vaf = vaf_list, dp = dp_list, gt = gt_list, ps = ps_list,
             caller = caller_name)
}

# Header-only provenance for the VCF(s) that supplied VAF/DP/GT/PS.
#
# Which file this is matters to a reader: when the pipeline ran several somatic callers
# through consensus, the FORMAT fields in the surviving VCF come from whichever caller won
# each merge, so they need not correspond to the `callers` column beside them. The report
# names the file rather than silently implying one caller — see the footnote in
# templates/sections/_smallvariants.qmd.
#
# `vcf_files` may be several paths (parse_caller_vcf() stacks them; ClairS splits matched
# mode into snvs.vcf.gz + indel.vcf.gz). Returns NULL when there is nothing to describe,
# matching the graceful-missing contract of the parsers above. Only the header is read —
# no record parsing — so this stays cheap on a multi-million-variant VCF.
vaf_provenance = function(vcf_files, sample_dir = NULL, sample_id = NULL) {
  if (is.null(vcf_files) || length(vcf_files) == 0) return(NULL)
  vcf_files = vcf_files[!is.na(vcf_files) & nzchar(vcf_files)]
  if (length(vcf_files) == 0) return(NULL)

  # Paths read better relative to the sample directory the caller passed in.
  rel = function(p) {
    if (is.null(sample_dir) || !nzchar(sample_dir)) return(p)
    root = sub("/+$", "", normalizePath(sample_dir, mustWork = FALSE))
    full = normalizePath(p, mustWork = FALSE)
    if (startsWith(full, paste0(root, "/"))) substring(full, nchar(root) + 2L) else p
  }

  # "##source=Clair3" and friends. Plenty of VCFs carry no source line at all, which is why
  # the footnote treats this as optional colour rather than the identifying fact.
  read_sources = function(p) {
    if (!file.exists(p)) return(character(0))
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    out = character(0)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) break
      if (!startsWith(line, "##")) break   # #CHROM or a malformed header ends the scan
      if (startsWith(line, "##source=")) out = c(out, sub("^##source=", "", line))
    }
    out
  }

  # Which sample column parse_caller_vcf() will have read. Worth naming in the footnote:
  # on a matched two-sample VCF the wrong pick reports the normal's VAF for every somatic
  # variant, and nothing else in the report would show it.
  read_sample = function(p) {
    if (!file.exists(p)) return(NA_character_)
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) return(NA_character_)
      if (startsWith(line, "#CHROM")) {
        header = strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(header) <= 9L) return(NA_character_)
        samples = header[10:length(header)]
        # Only worth reporting when there was a choice to get wrong.
        if (length(samples) == 1L) return(NA_character_)
        return(samples[suppressWarnings(pick_sample_column(samples, sample_id, p))])
      }
      if (!startsWith(line, "##")) return(NA_character_)
    }
  }
  chosen = unique(unlist(lapply(vcf_files, read_sample)))
  chosen = chosen[!is.na(chosen)]

  list(
    paths   = unname(vapply(vcf_files, rel, character(1))),
    sources = unique(unlist(lapply(vcf_files, read_sources))),
    sample  = chosen
  )
}

# Canonical variant key, used to join VEP annotation rows to the VCF they came from.
#
# VEP always reports an indel one base to the right of the VCF anchor, and writes the
# alleles either as the raw VCF pair or in its own trimmed form with a dash for the empty
# side — which of the two depends on the VEP version:
#
#   VCF record            VEP Uploaded_variation   allele notation
#   chr1 1871654 TG  T    chr1_1871655_TG/T        raw
#   chr1 14553006 G  GA   chr1_14553007_G/GA       raw
#   chr1 192936 GAATA G   chr1_192937_AATA/-       trimmed + dash
#   chr1 197487 A    AG   chr1_197488_-/G          trimmed + dash
#
# All four reconcile in one space: trimmed alleles (anchor base dropped, empty side
# written "-") at the VCF anchor position + 1. SNVs and equal-length MNVs have no anchor
# base and are keyed verbatim.
#
# This relies on `pos` coming from VEP's variant_id, which is consistently the shifted
# position — the Location column is not (see parse_vep_text()).
#
# `space` is "vcf" for records read from a VCF, "vep" for rows read from VEP output.
variant_key = function(chrom, pos, ref, alt, space = c("vcf", "vep")) {
  space = match.arg(space)
  ref = toupper(as.character(ref)); alt = toupper(as.character(alt))
  pos = as.integer(pos)

  trim = function(x) { t = substr(x, 2L, nchar(x)); ifelse(t == "", "-", t) }

  dash     = ref == "-" | alt == "-"   # already trimmed by VEP
  is_indel = dash | nchar(ref) != nchar(alt)

  # Raw allele pairs still need the anchor base dropped; dash forms are already trimmed.
  need_trim = is_indel & !dash

  # Only the VCF side needs shifting — VEP has already done it.
  key_pos = ifelse(space == "vcf" & is_indel, pos + 1L, pos)
  key_ref = ifelse(need_trim, trim(ref), ref)
  key_alt = ifelse(need_trim, trim(alt), alt)

  paste(chrom, key_pos, key_ref, key_alt, sep = "|")
}

# Classify SNV into 6 SBS mutation categories (C/T-ref normalised)
classify_mut = function(ref, alt) {
  comp = c(A = "T", T = "A", C = "G", G = "C")
  ref = toupper(ref); alt = toupper(alt)
  use_comp = !(ref %in% c("C", "T"))
  norm_ref = ifelse(use_comp, comp[ref], ref)
  norm_alt = ifelse(use_comp, comp[alt], alt)
  paste0(norm_ref, ">", norm_alt)
}

# Build the small-variant display table: canonical rows from the VEP annotation, with
# VAF / depth / phasing joined from the VCF that VEP annotated (see `somatic_vaf_vcf` in
# locate_outputs.R).
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all variants.
build_variant_table = function(vep_data, vaf_data, gene_panel = NULL) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)

  # Impact ranking for deduplication
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep_data[, impact_rank := impact_rank[impact]]
  vep_data[is.na(impact_rank), impact_rank := 5L]

  # Filter to gene panel (by gene symbol or Ensembl ID fallback)
  if (!is.null(gene_panel)) {
    if (length(gene_panel) > 0) {
      vep_data = vep_data[symbol %in% gene_panel | gene_id %in% gene_panel]
    } else {
      vep_data = vep_data[FALSE]  # Empty panel → empty result
    }
  }
  if (nrow(vep_data) == 0) return(data.table())

  # Keep best consequence per variant×gene (lowest impact rank)
  key_cols = c("chrom", "pos", "ref", "alt", "symbol")
  setorder(vep_data, impact_rank)
  vep_data = unique(vep_data, by = key_cols)

  # Join VEP rows to the VCF they were produced from, via the canonical key that
  # reconciles the two sides' indel representations (see variant_key()).
  #
  # The space is the parser's to declare, not ours to assume: parse_vep_text() returns
  # VEP-space coordinates (read from `variant_id`), parse_vep_vcf() returns the VCF's own
  # POS/REF/ALT. Hard-coding "vep" here shifted only one side of the key on the CSQ path,
  # so no indel could ever match and its VAF/DP/GT/PS came back silently NA.
  vep_space = if ("coord_space" %in% names(vep_data)) unique(vep_data$coord_space) else "vep"
  stopifnot(length(vep_space) == 1L, vep_space %in% c("vep", "vcf"))
  vep_data[, join_key := variant_key(chrom, pos, ref, alt, space = vep_space)]

  if (!is.null(vaf_data) && nrow(vaf_data) > 0) {
    vdt = vaf_data[, .(join_key = variant_key(chrom, pos, ref, alt, space = "vcf"),
                       vaf, dp, gt, ps)]
    vdt = unique(vdt, by = "join_key")
    vep_data = merge(vep_data, vdt, by = "join_key", all.x = TRUE)
  } else {
    vep_data[, `:=`(vaf = NA_real_, dp = NA_integer_,
                    gt  = NA_character_, ps = NA_character_)]
  }

  # Which caller reported each variant. A merged multi-caller VEP VCF states this outright
  # in INFO/CALLER; the VEP text format carries no per-variant caller, leaving this empty.
  if ("caller" %in% names(vep_data) && any(!is.na(vep_data$caller))) {
    vep_data[, callers := caller]
  } else {
    vep_data[, callers := ""]
  }

  # Mutation category for SNVs
  vep_data[nchar(ref) == 1 & nchar(alt) == 1,
           mut_cat := classify_mut(ref, alt)]

  display_cols = c("symbol", "chrom", "pos", "ref", "alt",
                   "consequence", "impact", "hgvsp",
                   "vaf", "dp", "gt", "ps",
                   "callers", "cosmic", "dbsnp", "sift", "polyphen")
  display_cols = display_cols[display_cols %in% names(vep_data)]
  vep_data[, ..display_cols]
}
