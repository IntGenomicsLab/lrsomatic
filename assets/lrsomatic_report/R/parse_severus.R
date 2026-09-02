suppressPackageStartupMessages({
  library(data.table)
})

# Severus writes both sides of a rearrangement as separate records (`_1` / `_2`,
# linked by INFO/MATE_ID) and tags single breakends "sBND", which fails a bare
# `== "BND"` test. Both are handled here, once, so that everything downstream
# — the SV table, the panel filter, the circos links and the SV count — sees one
# row per rearrangement carrying *both* breakend loci.
BND_SVTYPES = c("BND", "sBND")

# The BND-derived `svclass` values (see the fcase in parse_severus_somatic_records()).
# What these three share against DEL/DUP/INV/INS is what sv_display_columns() turns on:
# a junction is *two loci*, so it has no span and no size, while a contiguous type is
# *one span* whose two records are its own start and end. BND_CIRCOS_CLASSES
# (R/circos_bnd.R) is this set minus "single breakend", which has no partner locus and
# therefore no arc.
SV_JUNCTION_CLASSES = c("translocation", "intra-chr breakend", "single breakend")

# Panel-matching windows. Defined here as the single source of truth: they are
# used by sv_panel_hits() below and emitted to the report's client-side filter
# from templates/per_sample.qmd, so the R "Panel SVs" card and the JS row filter
# cannot drift apart.
SV_PANEL_WINDOW_BND   = 1e6   # distance from either breakend of a BND
SV_PANEL_WINDOW_OTHER = 1e5   # distance from the span of a DEL/DUP/INV/INS

# Read a VCF's data records, skipping the header. Returns an empty data.table
# rather than erroring for a header-only VCF.
.severus_read_vcf = function(vcf_file, col_names) {
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
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0) return(data.table())

  # Name the columns positionally rather than selecting a fixed count: a VCF with no
  # FORMAT/SAMPLE columns would otherwise fail the read and come back as zero SVs,
  # which reads identically to "this sample has none".
  n = min(ncol(dt), length(col_names))
  setnames(dt, seq_len(n), col_names[seq_len(n)])
  for (nm in setdiff(col_names, names(dt))) dt[, (nm) := NA_character_]
  dt[, ..col_names]
}

# One INFO sub-field per record, NA where the key is absent. Length-preserving:
# bare regmatches() silently drops non-matching elements, which is a length-mismatch
# error waiting for the first optional key (MATE_ID, END, SVLEN — i.e. most of them).
.info_val = function(info_vec, key) {
  m   = regexpr(paste0("(?:^|;)", key, "=([^;]+)"), info_vec, perl = TRUE)
  out = rep(NA_character_, length(info_vec))
  hit = m > 0
  if (any(hit)) out[hit] = sub(".*=", "", regmatches(info_vec, m))
  out
}

# Partner locus of a BND, from the ALT bracket notation — all four forms
# ("N[chr7:1234[", "]chr7:1234]N", "N]chr7:1234]", "[chr7:1234[N") and with or
# without a "chr" prefix on the contig, since a bare-contig reference would
# otherwise drop every breakend. A single breakend (sBND) has no bracket pair
# and yields NA.
.bnd_partner = function(alt) {
  pat = "^.*?[\\[\\]]([^\\[\\]:]+):([0-9]+)[\\[\\]].*$"
  has = grepl(pat, alt, perl = TRUE)
  list(
    chrom = ensure_chr_prefix(ifelse(has, sub(pat, "\\1", alt, perl = TRUE), NA_character_)),
    pos   = ifelse(has, suppressWarnings(as.integer(sub(pat, "\\2", alt, perl = TRUE))),
                   NA_integer_)
  )
}

# Parse the somatic Severus VCF into one row per rearrangement: both breakend loci,
# type, length and VAF. Mate records are collapsed (`_1` is side A, `_2` side B), so
# this is also what the SV count and the circos links are derived from.
#
# Returns: id, id_b, svtype, svclass, chrom_a, pos_a, chrom_b, pos_b, sv_len, vaf
# For non-BND types chrom_b/pos_b are the SV's own end; for a single breakend they
# are NA.
parse_severus_somatic_records = function(vcf_file) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(data.table())

  dt = .severus_read_vcf(vcf_file, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                     "FILTER", "INFO", "FORMAT", "SAMPLE1"))
  if (nrow(dt) == 0) return(data.table())

  dt[, CHROM := ensure_chr_prefix(CHROM)]
  dt[, SVTYPE := .info_val(INFO, "SVTYPE")]
  dt[, is_bnd := SVTYPE %in% BND_SVTYPES]

  # END and SVLEN are optional: a BND/INS-only VCF carries neither, so create the
  # columns unconditionally before anything reads them.
  dt[, `:=`(END = NA_integer_, SVLEN = NA_integer_)]
  dt[grepl("END=",   INFO, fixed = TRUE), END   := as.integer(.info_val(INFO, "END"))]
  dt[grepl("SVLEN=", INFO, fixed = TRUE), SVLEN := as.integer(.info_val(INFO, "SVLEN"))]
  dt[is.na(END), END := POS]

  # Partner locus. For a BND it comes from ALT, which names it even when the mate
  # record itself was filtered out of the VCF; for everything else it is the SV's
  # own end on the same contig.
  partner = .bnd_partner(dt$ALT)
  dt[, `:=`(chrom_b = fifelse(is_bnd, partner$chrom, CHROM),
            pos_b   = fifelse(is_bnd, partner$pos,   as.integer(END)))]

  dt[, mate_id := .info_val(INFO, "MATE_ID")]

  # One row per rearrangement. Group the two mate records on the unordered
  # {ID, MATE_ID} pair, falling back to the shared `_1`/`_2` stem for older VCFs
  # that carry no MATE_ID. Groups that are not exactly a pair (an unpaired `_1`,
  # an orphaned `_2`, an sBND) stay as singletons rather than being dropped.
  dt[, pair_key := fifelse(!is.na(mate_id),
                           paste(pmin(ID, mate_id), pmax(ID, mate_id), sep = "|"),
                           sub("_[12]$", "", ID))]
  dt[, side_rank := fifelse(grepl("_1$", ID), 1L, 2L)]
  dt[, n_in_pair := .N, by = pair_key]
  if (any(dt$n_in_pair > 2L)) {
    message("Severus VCF: ", sum(dt$n_in_pair > 2L),
            " records share a mate group with more than two members; kept unpaired.")
    dt[n_in_pair > 2L, pair_key := ID]
    dt[, n_in_pair := .N, by = pair_key]
  }
  setorder(dt, pair_key, side_rank, ID)
  dt = dt[dt[, .I[1L], by = pair_key]$V1]

  # The mate's record ID, which is what the per-side VEP annotation is keyed on.
  dt[, id_b := fifelse(n_in_pair == 2L & !is.na(mate_id), mate_id, NA_character_)]

  dt[, svclass := fcase(
    is_bnd & is.na(chrom_b),      "single breakend",
    is_bnd & chrom_b != CHROM,    "translocation",
    is_bnd,                       "intra-chr breakend",
    default = SVTYPE
  )]

  # VAF from FORMAT/SAMPLE1 (format string is uniform for Severus output, but split
  # format-group by format-group defensively, as in parse_caller_vcf())
  fmt_groups = unique(dt$FORMAT)
  vaf_list = rep(NA_real_, nrow(dt))
  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields = strsplit(fmt, ":", fixed = TRUE)[[1]]
    vaf_idx = match("VAF", fields)
    if (is.na(vaf_idx)) next
    split_s = strsplit(dt$SAMPLE1[idx_rows], ":", fixed = TRUE)
    vaf_list[idx_rows] = vapply(split_s, function(x)
      if (length(x) >= vaf_idx) suppressWarnings(as.numeric(x[vaf_idx])) else NA_real_,
      numeric(1))
  }
  dt[, VAF := vaf_list]

  dt[, .(id = ID, id_b, svtype = SVTYPE, svclass,
         chrom_a = CHROM, pos_a = POS, chrom_b, pos_b,
         sv_len = SVLEN, vaf = VAF)]
}

# Circos tracks from the collapsed records: BND links (one per rearrangement, so each
# arc is drawn once) and a non-BND track carrying colour/y-position.
severus_circos_tracks = function(records) {
  empty = list(translocations = data.table(), nontrans = data.table())
  if (is.null(records) || nrow(records) == 0) return(empty)

  is_bnd = records$svtype %in% BND_SVTYPES

  trans = records[is_bnd & !is.na(chrom_b) & !is.na(pos_b),
                  .(chrom = chrom_a, pos = pos_a, chrom2 = chrom_b, pos2 = pos_b)]
  # Records are already mate-collapsed; dedup on the unordered endpoint pair as well,
  # so a VCF without MATE_ID cannot draw the same arc twice.
  if (nrow(trans) > 0) {
    a = paste(trans$chrom, trans$pos, sep = ":")
    b = paste(trans$chrom2, trans$pos2, sep = ":")
    trans = unique(trans[, .link_key := paste(pmin(a, b), pmax(a, b))],
                   by = ".link_key")[, .link_key := NULL]
  }

  SV_COL  = c(INS = "#f97e02", DEL = "#020272", INV = "#e7cc02", DUP = "#e41a1c")
  SV_YPOS = c(INS = 1.0,       DEL = 0.66,      INV = 0.33,      DUP = 0.05)
  nontrans = records[!is_bnd,
    .(chrom = chrom_a, pos = pos_a, end = pos_b, svtype, svlen = sv_len,
      circos_pos = unname(SV_YPOS[svtype]), circos_col = unname(SV_COL[svtype]))]

  list(translocations = trans, nontrans = nontrans)
}

# Build the SV display table: mate-collapsed Severus records plus, as display context,
# the VEP gene symbol at each breakend (from the SV VEP VCF, `parse_vep_vcf()` in
# R/parse_smallvariants.R).
#
# The join is keyed on the VCF record ID, which is what links a VEP record back to the
# Severus record it was annotated from; a locus key is only used as a fallback for
# annotation files whose IDs don't overlap the Severus ones at all. Keying on the locus
# alone leaks annotations between distinct records that happen to start at the same base.
#
# These symbols are context, not the panel filter's key — see sv_panel_hits(). Whether
# VEP annotates a breakend at all is an invocation-dependent property of the sample
# (1.6%-90% of breakends across the samples measured), which is exactly why panel
# matching is done on coordinates instead.
build_sv_table_from_vep = function(records, vep_sv_vcf) {
  if (is.character(records)) records = parse_severus_somatic_records(records)
  if (is.null(records) || nrow(records) == 0) return(data.table())

  sv = copy(records)
  out_cols = c("id", "svclass", "svtype", "chrom_a", "pos_a", "gene_a",
               "chrom_b", "pos_b", "gene_b", "sv_len", "vaf",
               "consequence", "impact")

  vep = parse_vep_vcf(vep_sv_vcf)
  if (is.null(vep) || nrow(vep) == 0) {
    sv[, `:=`(gene_a = NA_character_, gene_b = NA_character_,
              consequence = NA_character_, impact = NA_character_)]
    return(sv[, ..out_cols])
  }

  # Highest-impact annotation first, so annot[1] below is the one that survives.
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep[, impact_rank := impact_rank[impact]]
  vep[is.na(impact_rank), impact_rank := 5L]
  setorder(vep, impact_rank)

  collapse_genes = function(x) {
    g = unique(x[!is.na(x) & nzchar(x)])
    if (length(g) == 0) NA_character_ else paste(g, collapse = ",")
  }

  by_id = if ("id" %in% names(vep) && any(!is.na(vep$id))) {
    vep[!is.na(id) & nzchar(id), .(gene = collapse_genes(symbol),
                                   consequence = consequence[1], impact = impact[1]),
        by = .(id)]
  } else NULL

  by_locus = vep[, .(gene = collapse_genes(symbol),
                     consequence = consequence[1], impact = impact[1]),
                 by = .(chrom, pos)]

  # Which key to use is decided once for the whole table, not per side: deciding it per
  # side would silently locus-match side B (where BND mate IDs are the only IDs) while
  # ID-matching side A.
  use_id = !is.null(by_id) && any(c(sv$id, sv$id_b) %in% by_id$id)

  # Severus records that the annotation side filtered out (parse_vep_vcf() keeps only
  # FILTER PASS/.) simply get no symbols — that is not a join bug.
  side_annot = function(ids, chroms, positions) {
    src = if (use_id) by_id else by_locus
    m   = if (use_id) match(ids, by_id$id)
          else match(paste(chroms, positions), paste(by_locus$chrom, by_locus$pos))
    res = data.table(gene = NA_character_, consequence = NA_character_,
                     impact = NA_character_)[rep(1L, length(chroms))]
    hit = !is.na(m)
    if (any(hit)) res[hit, `:=`(gene = src$gene[m[hit]],
                                consequence = src$consequence[m[hit]],
                                impact = src$impact[m[hit]])]
    res
  }

  a = side_annot(sv$id,   sv$chrom_a, sv$pos_a)
  b = side_annot(sv$id_b, sv$chrom_b, sv$pos_b)

  sv[, `:=`(gene_a = a$gene, gene_b = b$gene)]
  # Consequence/impact describe the rearrangement, taken from the higher-impact side.
  rank_of = function(x) { r = unname(impact_rank[x]); r[is.na(r)] = 5L; r }
  use_b   = rank_of(b$impact) < rank_of(a$impact)
  sv[, `:=`(consequence = fifelse(use_b, b$consequence, a$consequence),
            impact      = fifelse(use_b, b$impact,      a$impact))]

  sv[, ..out_cols]
}

# Human-readable location, size and gene columns for the SV table.
#
# The stored schema is bedpe-shaped (chrom_a/pos_a + chrom_b/pos_b) because that is the
# only shape a translocation fits, and everything downstream reads it — sv_panel_hits(),
# the client-side filter, bnd_links(). But it is the wrong *reading* for most rows: a
# DEL's two records are its own start and end, not two partners, and its gene_a/gene_b
# are just the genes at either edge of one span. So the raw columns stay (hidden in the
# rendered table, and still exported) and these are shown instead.
#
# The rule is `svclass`, not "same chromosome": parse_severus_somatic_records() sets
# chrom_b = chrom_a for DEL/DUP/INV/INS, so chrom_b == chrom_a is equally true of an
# intra-chromosomal breakend — and there the two loci are a junction, not a span the
# variant swallows. Keying off the contig would hand an intra-chr breakend a start, an
# end and a size it does not have.
#
# Returns a data.table of five columns, one row per row of `sv_table`, same order:
#   locus      the display string (see the fcase below)
#   locus_sort sort key for `locus`, which as a display string sorts lexically
#   size       fmt_bp() of size_bp; "" for a junction
#   size_bp    the number behind `size`, and its sort key
#   genes      merged for a span, side-labelled for a junction
# `chrom_levels` orders the sort key: pass the report's plotted chromosome vector from
# chromosomes_for_sex(), which is already in natural order. Contigs outside it sort last.
sv_display_columns = function(sv_table, chrom_levels = NULL) {
  empty = data.table(locus = character(), locus_sort = character(),
                     size = character(), size_bp = numeric(), genes = character())
  if (is.null(sv_table) || nrow(sv_table) == 0) return(empty)

  d = as.data.table(sv_table)
  need = c("svclass", "chrom_a", "pos_a", "chrom_b", "pos_b", "sv_len",
           "gene_a", "gene_b")
  absent = setdiff(need, names(d))
  if (length(absent) > 0)
    stop("sv_display_columns(): sv_table is missing ", paste(absent, collapse = ", "))

  is_junction = d$svclass %in% SV_JUNCTION_CLASSES
  pa = suppressWarnings(as.numeric(d$pos_a))
  pb = suppressWarnings(as.numeric(d$pos_b))

  # formatC(), not format(): format() is vectorised to a *common* width and would pad
  # every position out to the longest one in the table.
  bp = function(x) formatC(x, format = "d", big.mark = ",")
  at = function(chrom, pos) paste0(chrom, ":", bp(pos))

  locus = fcase(
    is.na(pa) | is.na(d$chrom_a),   NA_character_,
    d$svclass == "single breakend", paste0(at(d$chrom_a, pa), " (unpaired)"),
    is_junction,
      # An arrow, not a dash: these two loci are joined, they do not bound a span. Both
      # sides are named even for an intra-chromosomal junction, so a per-column search
      # for a chromosome matches it on either side.
      fifelse(is.na(pb) | is.na(d$chrom_b),
              paste0(at(d$chrom_a, pa), " (unpaired)"),
              paste0(at(d$chrom_a, pa),
                     fifelse(d$svclass == "translocation", " → ", " ↔ "),
                     at(d$chrom_b, pb))),
    # A contiguous type: one span. An INS has END == POS, so it is a single point.
    default = fifelse(is.na(pb) | pb <= pa,
                      at(d$chrom_a, pa),
                      paste0(at(d$chrom_a, pa), "–", bp(pb)))
  )

  # SVLEN where Severus wrote one — an INS's length is *not* its span — else the span.
  len = suppressWarnings(as.numeric(d$sv_len))
  size_bp = fifelse(is_junction, NA_real_, fifelse(!is.na(len), abs(len), pb - pa))
  # as.character() is load-bearing: fmt_bp() is built on ifelse(), which returns a vector
  # typed after its *test*, so an all-NA size_bp — every row a junction, i.e. a BND-only
  # SV table — comes back logical and the fifelse() below would reject the type mismatch.
  size = as.character(fmt_bp(size_bp))

  split_genes = function(x) {
    if (is.na(x) || !nzchar(x)) return(character())
    g = trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    unique(g[nzchar(g)])
  }
  # A middot rather than the whitespace a mock-up would use: DT escapes cell content, so
  # &nbsp; is unavailable, and HTML collapses a run of spaces to one.
  genes = vapply(seq_len(nrow(d)), function(i) {
    a = split_genes(d$gene_a[i]); b = split_genes(d$gene_b[i])
    if (!is_junction[i]) return(paste(unique(c(a, b)), collapse = ", "))
    sides = c(if (length(a)) paste0("A: ", paste(a, collapse = ",")),
              if (length(b)) paste0("B: ", paste(b, collapse = ",")))
    paste(sides, collapse = " · ")
  }, character(1))

  rank = if (length(chrom_levels)) match(d$chrom_a, chrom_levels)
         else rep(NA_integer_, nrow(d))
  rank[is.na(rank)] = 999L
  locus_sort = fifelse(is.na(pa), "", sprintf("%03d:%011.0f", rank, pa))

  data.table(locus      = fifelse(is.na(locus), "", locus),
             locus_sort = locus_sort,
             size       = fifelse(is.na(size), "", size),
             size_bp    = size_bp,
             genes      = genes)
}

# Which panel genes each SV hits, and how it hit them.
#
# Coordinate-carrying panels (see load_gene_panel()) are matched on position: a panel
# interval within `bnd_window` of either breakend of a BND, or within `other_window` of
# the span of any other type. That needs no annotation on the SV row, which is the point
# — VEP breakend coverage is sample-dependent, so symbol matching hides breakends on
# samples whose VEP run didn't annotate them.
#
# A symbol-only panel falls back to a direct hit on either side's VEP symbol: today's
# behaviour, but reading both breakends instead of one.
#
# Returns one character label per row of `sv_table` ("" = no hit). Each hit reads
# "GENE (side, how)": `side` is the breakend it matched on (A/B, or "span" for a
# contiguous type), and `how` is "direct" when the locus overlaps the gene interval, or
# the distance to it otherwise. Without that second token a breakend most of a megabase
# away reads exactly like one inside the gene — the BND window is wide, and proximity is
# not disruption. The "Panel SVs" summary card counts the non-empty entries; the
# client-side filter applies the same test with the same windows and builds the same
# labels.
#
# Several panels can be active at once, so `panels` is either a single panel object or a
# named list of them; a hit against any of them counts (union). When more than one is
# active each label gains a " [name]" suffix saying which panel matched — with one panel
# the labels are exactly what they always were, since there is nothing to disambiguate.
# svPanelHits() in the panel-js-data chunk of templates/per_sample.qmd mirrors this rule
# and must change with it.
sv_panel_hits = function(sv_table, panels,
                         bnd_window = SV_PANEL_WINDOW_BND,
                         other_window = SV_PANEL_WINDOW_OTHER) {
  n = if (is.null(sv_table)) 0L else nrow(sv_table)
  if (n == 0) return(character(0))

  ps = .as_panel_list(panels)
  if (length(ps) == 0) return(rep("", n))
  if (length(ps) == 1) return(.sv_panel_hits_one(sv_table, ps[[1]], "", bnd_window, other_window))

  per_panel = lapply(names(ps), function(nm)
    .sv_panel_hits_one(sv_table, ps[[nm]], paste0(" [", nm, "]"), bnd_window, other_window))
  vapply(seq_len(n), function(i) {
    parts = unique(unlist(lapply(per_panel, `[[`, i), use.names = FALSE))
    paste(parts[nzchar(parts)], collapse = ", ")
  }, character(1))
}

# A panel object is a plain list, and so is a list of panels — tell them apart by the
# fields load_gene_panel() always sets rather than by class.
.as_panel_list = function(panels) {
  if (is.null(panels) || length(panels) == 0) return(list())
  if (!is.null(panels$genes) || !is.null(panels$has_coords)) {
    nm = if (!is.null(panels$name)) as.character(panels$name)[1] else "panel"
    return(setNames(list(panels), nm))
  }
  ps = panels[!vapply(panels, is.null, logical(1))]
  if (is.null(names(ps))) names(ps) = paste0("panel", seq_along(ps))
  ps
}

# One panel's per-row labels. `tag` is appended to every label ("" for a lone panel).
.sv_panel_hits_one = function(sv_table, panel, tag = "",
                              bnd_window = SV_PANEL_WINDOW_BND,
                              other_window = SV_PANEL_WINDOW_OTHER) {
  n = nrow(sv_table)
  if (is.null(panel)) return(rep("", n))

  label = function(genes, side, how) paste0(genes, " (", side, ", ", how, ")", tag)

  if (!isTRUE(panel$has_coords)) {
    symbols = toupper(panel$genes)
    hit_side = function(col) {
      if (!col %in% names(sv_table)) return(rep("", n))
      vapply(sv_table[[col]], function(cell) {
        if (is.na(cell) || !nzchar(cell)) return("")
        g = trimws(unlist(strsplit(toupper(as.character(cell)), "[;,]+")))
        g = unique(g[g %in% symbols])
        if (length(g) == 0) "" else paste(g, collapse = ",")
      }, character(1), USE.NAMES = FALSE)
    }
    ha = hit_side("gene_a"); hb = hit_side("gene_b")
    # A VEP symbol sits on the breakend itself, so a symbol hit is always direct — and
    # there are no coordinates on this path to measure anything else from.
    return(vapply(seq_len(n), function(i) {
      parts = c(if (nzchar(ha[i])) label(ha[i], "A", "direct"),
                if (nzchar(hb[i])) label(hb[i], "B", "direct"))
      paste(parts, collapse = ", ")
    }, character(1)))
  }

  iv = panel_intervals(panel)
  if (is.null(iv) || nrow(iv) == 0) return(rep("", n))

  is_bnd = sv_table$svtype %in% BND_SVTYPES
  # `start`/`end` are padded by the window and are what foverlaps() matches on. `qlo`/`qhi`
  # carry the *unpadded* locus through the join: the gap to a gene has to be measured from
  # where the breakend actually is, not from the edge of the window around it.
  q = rbindlist(list(
    # Each breakend of a BND gets its own window; the two sides can be on
    # different contigs, so they cannot be one interval.
    data.table(row = which(is_bnd), side = "A",
               chrom = sv_table$chrom_a[is_bnd],
               start = sv_table$pos_a[is_bnd] - bnd_window,
               end   = sv_table$pos_a[is_bnd] + bnd_window,
               qlo   = sv_table$pos_a[is_bnd],
               qhi   = sv_table$pos_a[is_bnd]),
    data.table(row = which(is_bnd), side = "B",
               chrom = sv_table$chrom_b[is_bnd],
               start = sv_table$pos_b[is_bnd] - bnd_window,
               end   = sv_table$pos_b[is_bnd] + bnd_window,
               qlo   = sv_table$pos_b[is_bnd],
               qhi   = sv_table$pos_b[is_bnd]),
    # Other types are contiguous: one window around the whole span.
    data.table(row = which(!is_bnd), side = "span",
               chrom = sv_table$chrom_a[!is_bnd],
               start = pmin(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]) - other_window,
               end   = pmax(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]) + other_window,
               qlo   = pmin(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]),
               qhi   = pmax(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]))
  ))
  q = q[!is.na(chrom) & !is.na(start) & !is.na(end)]
  if (nrow(q) == 0) return(rep("", n))
  q[, start := pmax(as.numeric(start), 0)]
  q[, end   := as.numeric(end)]
  q[, `:=`(qlo = as.numeric(qlo), qhi = as.numeric(qhi))]
  iv = copy(iv)
  iv[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  setkey(iv, chrom, start, end)

  ov = data.table::foverlaps(q, iv, by.x = c("chrom", "start", "end"),
                             type = "any", nomatch = NULL)
  if (nrow(ov) == 0) return(rep("", n))

  # `start`/`end` are the gene's own interval here — foverlaps() renamed the padded query
  # window to i.start/i.end. A zero gap means the locus lands inside the gene.
  ov[, gap := pmax(0, pmax(start - qhi, qlo - end))]
  ov[, how := fifelse(gap == 0, "direct", fmt_bp(gap))]

  per_row = ov[, .(hit = paste(unique(label(gene, side, how)), collapse = ", ")), by = row]
  out = rep("", n)
  out[per_row$row] = per_row$hit
  out
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

# Build the SV display table from the gene-annotated Severus TSV — the fallback for
# pipelines that produce it instead of a VEP SV VCF. One row per SV, mapped onto the
# same column contract as the VEP path, so the panel filter, the counts and the section
# template have one schema to know about (and so this path gets coordinate matching too).
build_sv_table = function(sv_tsv) {
  if (is.null(sv_tsv) || nrow(sv_tsv) == 0) return(data.table())

  col = function(nm, default = NA) if (nm %in% names(sv_tsv)) sv_tsv[[nm]] else default
  chrom_a = ensure_chr_prefix(as.character(col("START_CHROM", NA_character_)))
  chrom_b = ensure_chr_prefix(as.character(col("END_CHROM",   NA_character_)))
  svtype  = as.character(col("SVTYPE", NA_character_))

  data.table(
    id      = as.character(col("ID", NA_character_)),
    svclass = fcase(
      is.na(chrom_a) | is.na(chrom_b), svtype,
      svtype %in% BND_SVTYPES & chrom_a != chrom_b, "translocation",
      svtype %in% BND_SVTYPES, "intra-chr breakend",
      default = svtype
    ),
    svtype  = svtype,
    chrom_a = chrom_a,
    pos_a   = suppressWarnings(as.integer(col("START_POS", NA_integer_))),
    # The TSV's gene hits are span-level, with no per-breakend split — hence side A only.
    gene_a  = as.character(if ("NHL_GENE_HITS" %in% names(sv_tsv)) sv_tsv$NHL_GENE_HITS
                           else col("COSMIC_GENE_HITS", NA_character_)),
    chrom_b = chrom_b,
    pos_b   = suppressWarnings(as.integer(col("END_POS", NA_integer_))),
    gene_b  = NA_character_,
    sv_len  = suppressWarnings(as.integer(col("SV_LEN", NA_integer_))),
    vaf     = suppressWarnings(as.numeric(col("VAF", NA_real_))),
    consequence = as.character(col("DETAILED_TYPE", NA_character_)),
    impact  = NA_character_
  )
}
