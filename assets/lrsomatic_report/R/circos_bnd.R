suppressPackageStartupMessages({
  library(data.table)
})

# Breakend circos: a second, BND-only plot drawn over just the chromosomes a breakend
# touches, cross-linked to the SV table so selecting a row highlights its arc.
#
# R selects the data; the browser draws it. That split exists because the plot has to
# re-lay-out when the gene panel changes: which chromosomes get a sector, and therefore
# every sector's angular width, is a function of the *filtered* link set. A server-drawn
# SVG bakes those angles in, and its arcs are flattened point lists carrying no genomic
# coordinates, so no client script can recompute them — the previous version could only
# dim arcs in place, leaving every chromosome on the plot whatever the filter said.
#
# This also retires the sentinel-colour tagging that the inlined-SVG design needed
# (svglite emits no ids, so each object had to be drawn in a unique colour and the markup
# rewritten afterwards), along with its "a mistagged arc disables highlighting" fallback:
# nodes the client creates itself carry their own ids by construction.
#
# See assets/js/bnd_circos.js for the drawing half.

# Above this the payload stops being worth its bytes: real samples on this cluster carry
# 50-115 arcs after mate collapse, but an unfiltered run can carry ~13k.
BND_CIRCOS_MAX_LINKS = 2000L

# The BND classes with two loci to draw an arc between; svclass is set by
# parse_severus_somatic_records() / build_sv_table(). A "single breakend" has no partner
# locus and so cannot be drawn at all.
BND_CIRCOS_CLASSES = c("translocation", "intra-chr breakend")

# ---- Data selection ------------------------------------------------------

# The rearrangements the circos can draw: BNDs with both loci known, on chromosomes the
# report is plotting. Callers pass the result to both bnd_panel_genes() and
# bnd_circos_data(), so the gene track cannot annotate an arc that was never drawn.
bnd_links = function(sv_table, chromosomes) {
  cols = c("id", "svclass", "chrom_a", "pos_a", "chrom_b", "pos_b")
  empty = data.table(id = character(), svclass = character(),
                     chrom_a = character(), pos_a = integer(),
                     chrom_b = character(), pos_b = integer())
  if (is.null(sv_table) || nrow(sv_table) == 0 || !all(cols %in% names(sv_table)))
    return(empty)

  b = as.data.table(sv_table)[, ..cols]
  b = b[svclass %in% BND_CIRCOS_CLASSES &
        !is.na(chrom_a) & !is.na(pos_a) & !is.na(chrom_b) & !is.na(pos_b) &
        chrom_a %in% chromosomes & chrom_b %in% chromosomes]
  if (nrow(b) == 0) return(empty)

  # Deterministic order, so the payload is byte-identical across renders of the same data.
  setorder(b, chrom_a, pos_a, chrom_b, pos_b, id)
  b[]
}

# The panel genes to label: for every coordinate-carrying panel, those whose interval
# falls within `window` of either breakend of a drawn arc.
#
# This is the same test, with the same window, that sv_panel_hits() applies to fill the
# table's panel_hit column — deliberately, because the report only shows a gene here once
# it appears in that column. Reimplementing the predicate would let the two drift.
#
# A symbol-only panel (and the custom paste-in panel) carries no coordinates, so it can
# produce panel_hit labels but no gene bodies.
bnd_panel_genes = function(links, all_panels, window = SV_PANEL_WINDOW_BND) {
  empty = data.table(chrom = character(), start = integer(), end = integer(),
                     gene = character(), panels = character())
  if (is.null(links) || nrow(links) == 0) return(empty)
  if (is.null(all_panels) || length(all_panels) == 0) return(empty)

  q = rbindlist(list(
    data.table(chrom = links$chrom_a, pos = as.numeric(links$pos_a)),
    data.table(chrom = links$chrom_b, pos = as.numeric(links$pos_b))
  ))
  q = unique(q[!is.na(chrom) & !is.na(pos)])
  if (nrow(q) == 0) return(empty)
  q[, `:=`(start = pmax(pos - window, 0), end = pos + window)]

  out = list()
  for (nm in names(all_panels)) {
    iv = panel_intervals(all_panels[[nm]])
    if (is.null(iv) || nrow(iv) == 0) next
    iv = copy(iv)[, `:=`(start = as.numeric(start), end = as.numeric(end))]
    setkey(iv, chrom, start, end)
    ov = data.table::foverlaps(q[, .(chrom, start, end)], iv,
                               by.x = c("chrom", "start", "end"),
                               type = "any", nomatch = NULL)
    if (nrow(ov) == 0) next
    out[[nm]] = unique(ov[, .(chrom, start, end, gene, panel = nm)])
  }
  if (length(out) == 0) return(empty)

  g = rbindlist(out)
  g = g[, .(panels = paste(sort(unique(panel)), collapse = ",")),
        by = .(chrom, start, end, gene)]
  setorder(g, chrom, start, gene)
  g[, `:=`(start = as.integer(start), end = as.integer(end))]
  g[, .(chrom, start, end, gene, panels)]
}


# ---- Payload for the client-side plot -------------------------------------

# Everything the browser needs to draw the breakend circos at any filter state.
#
# Returns list(data, n_links, n_genes, chroms, reason). `data` is NULL when there is
# nothing to draw or too much of it, with `reason` saying which; the template shows that
# through section_notice() instead of an empty ring.
#
# `chroms` here is the *unfiltered* superset — every chromosome some drawable breakend
# touches. The client narrows it per filter; R only needs to bound the payload.
bnd_circos_data = function(links, genes = NULL, cytobands, chrom_lengths, chromosomes) {

  fail = function(reason) list(data = NULL, n_links = 0L, n_genes = 0L,
                               chroms = character(0), reason = reason)

  if (is.null(links) || nrow(links) == 0)
    return(fail("No breakends with two mapped loci — nothing to draw arcs between."))
  if (nrow(links) > BND_CIRCOS_MAX_LINKS)
    return(fail(sprintf(paste("%s breakend arcs is more than this plot can show",
                              "(limit %s). Select a gene panel above for a filtered",
                              "report, or read the table below."),
                        nrow(links), BND_CIRCOS_MAX_LINKS)))

  b = copy(as.data.table(links))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]
  chroms = chroms[chroms %in% unique(cytobands$chrom)]
  b = b[chrom_a %in% chroms & chrom_b %in% chroms]
  if (nrow(b) == 0 || length(chroms) == 0)
    return(fail("No breakends on the chromosomes this report plots."))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]

  lens = chrom_lengths[chroms]
  # A locus past the end of its sector would be drawn outside the ring; clamp, as the
  # circlize version did, rather than dropping the row and losing its table cross-link.
  b[, `:=`(pos_a = pmin(pmax(as.numeric(pos_a), 1), lens[chrom_a]),
           pos_b = pmin(pmax(as.numeric(pos_b), 1), lens[chrom_b]))]

  cy = as.data.table(cytobands)[chrom %in% chroms]
  setorder(cy, chrom, start)

  g = if (is.null(genes)) data.table() else as.data.table(genes)[chrom %in% chroms]
  if (nrow(g) > 0) setorder(g, chrom, start, gene)

  list(
    data = list(
      chromosomes = js_vec(chroms),
      lengths     = js_vec(as.numeric(lens)),
      cytobands   = js_rows(cy, c("chrom", "start", "end", "stain")),
      links       = js_rows(b,  c("id", "svclass", "chrom_a", "pos_a", "chrom_b", "pos_b")),
      genes       = if (nrow(g) > 0) js_rows(g, c("chrom", "start", "end", "gene", "panels"))
                    else "[]"
    ),
    n_links = nrow(b),
    n_genes = nrow(g),
    chroms  = chroms,
    reason  = NULL
  )
}

# The <script> that hands the payload to the client, as one window.BND_DATA literal.
bnd_circos_script = function(res) {
  if (is.null(res$data)) return("")
  d = res$data
  paste0(
    "<script>window.BND_DATA={",
    "chromosomes:", d$chromosomes, ",",
    "lengths:",     d$lengths,     ",",
    "cytobands:",   d$cytobands,   ",",
    "links:",       d$links,       ",",
    "genes:",       d$genes,
    "};</script>")
}
