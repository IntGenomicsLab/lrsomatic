suppressPackageStartupMessages({
  library(data.table)
})

# Type-stable: fread() reads a bare-contig CHROM column ("1", "2", ...) as integer, and
# ifelse() on an all-NA test returns logical NA — both of which break every downstream
# string comparison, so coerce on the way in and on the way out.
ensure_chr_prefix = function(x) {
  x   = as.character(x)
  out = as.character(ifelse(startsWith(x, "chr"), x, paste0("chr", x)))
  out[is.na(x)] = NA_character_
  out
}

strip_chr_prefix = function(x) {
  sub("^chr", "", x)
}

# Parse VEP "Extra" key=value semicolon-delimited field into a named character vector
parse_extra_kv = function(extra_string) {
  if (is.na(extra_string) || extra_string == "" || extra_string == "-") return(character(0))
  pairs = strsplit(extra_string, ";", fixed = TRUE)[[1]]
  kv = strsplit(pairs, "=", fixed = TRUE)
  keys = vapply(kv, `[`, character(1), 1)
  vals = vapply(kv, function(x) if (length(x) >= 2) paste(x[-1], collapse = "=") else "", character(1))
  setNames(vals, keys)
}

# Vectorised: extract one key from VEP Extra column for each row
extract_extra_key = function(extra_vec, key) {
  vapply(extra_vec, function(x) {
    kv = parse_extra_kv(x)
    if (key %in% names(kv)) kv[[key]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# ---- Gene panels ---------------------------------------------------------
#
# A panel is a plain list — not a data.table — because it round-trips through
# Quarto's execute_params as YAML (see bin/render_report.R):
#
#   list(name, path, reference, has_coords, genes,
#        chrom, start, end, interval_gene)   # the last four only when has_coords
#
# `genes` is deduplicated, for symbol matching; the interval vectors are parallel and
# not deduplicated, because one symbol can legitimately carry several loci.
#
# `has_coords` selects the SV matching mode: a panel carrying chrom/start/end is
# matched against SV breakend coordinates with a window (sv_panel_hits()), while a
# symbol-only panel can only be matched against the per-side VEP symbols. Small
# variants always match on `genes` — VEP annotates SNVs reliably, and per-gene
# symbols are the right semantics there — so one file serves both tables.

# Canonical reference names. Panel coordinates are only valid for one reference, so
# a coordinate-carrying panel has to declare which, and the declaration is compared
# against the reference the report was rendered for.
normalise_reference_name = function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  x = tolower(trimws(x))
  if (x %in% c("t2t", "chm13", "chm13v2", "chm13v2.0", "t2t-chm13")) return("t2t")
  if (x %in% c("hg38", "grch38", "hg38-noalt")) return("hg38")
  x
}

# Reference suffix in a builtin panel filename ("lymphoid.hg38.tsv" -> "hg38").
# Returns list(name, reference); reference is NA for an unsuffixed file.
.split_panel_filename = function(path) {
  stem  = tools::file_path_sans_ext(basename(path))
  parts = strsplit(stem, ".", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    ref = normalise_reference_name(parts[length(parts)])
    if (!is.na(ref) && ref %in% c("t2t", "hg38"))
      return(list(name = paste(parts[-length(parts)], collapse = "."), reference = ref))
  }
  list(name = stem, reference = NA_character_)
}

# Read the leading "#"-comment block of a panel TSV, which may declare the reference
# the coordinates belong to ("# reference: hg38").
.panel_header = function(path) {
  lines = readLines(path, warn = FALSE)
  n_comment = 0L
  declared  = NA_character_
  for (ln in lines) {
    if (!startsWith(ln, "#")) break
    n_comment = n_comment + 1L
    m = regmatches(ln, regexpr("^#\\s*reference\\s*:\\s*\\S+", ln))
    if (length(m) > 0) declared = sub("^#\\s*reference\\s*:\\s*", "", m)
  }
  list(n_comment = n_comment, reference = declared)
}

# Load a gene panel TSV. Required: a `gene` column (or a single unnamed column of
# symbols). Optional but all-or-nothing: `chrom` (or `chr`), `start`, `end` — a file
# carrying some but not all three is malformed and errors rather than silently
# downgrading to symbol matching.
#
# `reference` is the reference the report is being rendered for. A coordinate-carrying
# panel that declares a different one errors; one that declares none loads with
# reference "" ("unverified" in the section footnote) rather than being guessed at.
load_gene_panel = function(path, reference = NULL) {
  if (!file.exists(path)) stop("Gene panel file not found: ", path)

  hdr = .panel_header(path)
  dt = tryCatch(
    fread(path, header = TRUE, sep = "\t", fill = TRUE, skip = hdr$n_comment),
    error = function(e) fread(path, header = FALSE, sep = "\t", fill = TRUE,
                              skip = hdr$n_comment)
  )
  if (nrow(dt) == 0 && ncol(dt) == 0) stop("Gene panel file is empty: ", path)

  # copy(): setnames() rewrites the names vector in place, which would otherwise
  # clobber this reference to it.
  orig_names = copy(names(dt))
  setnames(dt, tolower(names(dt)))
  setnames(dt, old = c("chr", "chromosome"), new = c("chrom", "chrom"), skip_absent = TRUE)

  if ("gene" %in% names(dt)) {
    genes = as.character(dt[["gene"]])
  } else if (ncol(dt) == 1) {
    # A headerless one-column list: fread consumed the first symbol as the column
    # name, so put it back — with its original casing — rather than dropping it.
    genes = c(orig_names[1], as.character(dt[[1]]))
  } else {
    genes = as.character(dt[[1]])
  }
  genes = trimws(genes)
  keep  = nzchar(genes) & !is.na(genes) & genes != "-"

  coord_cols = c("chrom", "start", "end")
  present    = intersect(coord_cols, names(dt))
  if (length(present) > 0 && length(present) < 3) {
    stop("Gene panel ", path, " carries coordinate column(s) ",
         paste(present, collapse = ", "), " but not all of ",
         paste(coord_cols, collapse = ", "),
         ". Supply all three, or none for symbol-only matching.")
  }
  has_coords = length(present) == 3

  # A `reference` column is an alternative to the "# reference:" comment line.
  declared = hdr$reference
  if (is.na(declared) && "reference" %in% names(dt)) {
    vals = unique(trimws(as.character(dt[["reference"]])))
    vals = vals[nzchar(vals) & !is.na(vals)]
    if (length(vals) > 1)
      stop("Gene panel ", path, " declares more than one reference: ",
           paste(vals, collapse = ", "))
    if (length(vals) == 1) declared = vals
  }
  declared_norm = normalise_reference_name(declared)
  want          = normalise_reference_name(reference)

  # Only coordinate panels are reference-specific; a symbol-only panel is agnostic.
  if (has_coords && !is.na(declared_norm) && !is.na(want) && declared_norm != want) {
    stop("Gene panel ", path, " declares reference '", declared_norm,
         "' but the report is being rendered against '", want,
         "'. Panel coordinates are only valid for the reference they were built on.")
  }

  out = list(
    name       = .split_panel_filename(path)$name,
    path       = path,
    reference  = if (is.na(declared_norm)) "" else declared_norm,
    has_coords = has_coords,
    genes      = unique(genes[keep])
  )

  if (has_coords) {
    chrom = ensure_chr_prefix(trimws(as.character(dt[["chrom"]])))
    start = suppressWarnings(as.integer(dt[["start"]]))
    end   = suppressWarnings(as.integer(dt[["end"]]))
    bad   = keep & (is.na(chrom) | !nzchar(chrom) | is.na(start) | is.na(end))
    if (any(bad))
      stop("Gene panel ", path, " has missing or non-numeric coordinates for: ",
           paste(utils::head(genes[bad], 5), collapse = ", "),
           if (sum(bad) > 5) paste0(" (and ", sum(bad) - 5L, " more)") else "")
    out$chrom = chrom[keep]
    out$start = start[keep]
    out$end   = end[keep]
    # `genes` is deduplicated for symbol matching; the interval vectors are not,
    # because one symbol can legitimately carry several loci.
    out$interval_gene = genes[keep]
  }

  out
}

# Panel intervals as a keyed data.table, or NULL for a symbol-only panel.
panel_intervals = function(panel) {
  if (is.null(panel) || !isTRUE(panel$has_coords)) return(NULL)
  dt = data.table(gene  = as.character(panel$interval_gene),
                  chrom = as.character(panel$chrom),
                  start = as.integer(panel$start),
                  end   = as.integer(panel$end))
  setkey(dt, chrom, start, end)
  dt
}

# Is a --gene-panel argument the "no filtering" sentinel?
is_no_gene_panel = function(panel_arg) {
  is.null(panel_arg) || length(panel_arg) != 1 || is.na(panel_arg) ||
    identical(tolower(trimws(panel_arg)), "none")
}

# Path of a builtin panel, preferring the variant built for `reference`
# ("lymphoid.hg38.tsv") over an unsuffixed one ("lymphoid.tsv"). NULL if neither exists.
builtin_panel_path = function(assets_dir, name, reference = NULL) {
  ref = normalise_reference_name(reference)
  dir = file.path(assets_dir, "gene_lists")
  candidates = c(if (!is.na(ref)) file.path(dir, paste0(name, ".", ref, ".tsv")),
                 file.path(dir, paste0(name, ".tsv")))
  hit = candidates[file.exists(candidates)]
  if (length(hit) > 0) hit[1] else NULL
}

# Resolve a --gene-panel arg: the "none" sentinel (no filtering, returns NULL),
# a builtin name ("lymphoid"), or a path to a TSV. A value that is neither is an
# error rather than a silent fall-through to unfiltered output.
resolve_gene_panel = function(panel_arg, assets_dir, reference = NULL) {
  if (is_no_gene_panel(panel_arg)) return(NULL)
  builtin = builtin_panel_path(assets_dir, panel_arg, reference)
  if (!is.null(builtin)) return(load_gene_panel(builtin, reference))
  if (file.exists(panel_arg)) return(load_gene_panel(panel_arg, reference))
  stop("Gene panel not found (tried builtin '", panel_arg, "' and as file path)")
}

# Load all gene panels from assets/gene_lists/*.tsv, resolving reference-specific
# files ("lymphoid.hg38.tsv" / "lymphoid.t2t.tsv") to one selectable "lymphoid" entry.
# A panel that only ships for other references is skipped — offering it would mean
# matching coordinates from the wrong genome.
# Returns a named list of panel objects (see load_gene_panel()).
load_all_gene_panels = function(assets_dir, reference = NULL) {
  tsv_files = Sys.glob(file.path(assets_dir, "gene_lists", "*.tsv"))
  if (length(tsv_files) == 0) return(list())

  meta = lapply(tsv_files, .split_panel_filename)
  ref  = normalise_reference_name(reference)

  panels = list()
  for (nm in unique(vapply(meta, `[[`, character(1), "name"))) {
    idx  = which(vapply(meta, `[[`, character(1), "name") == nm)
    refs = vapply(meta[idx], function(m) m$reference, character(1))
    pick = if (!is.na(ref) && any(refs == ref, na.rm = TRUE)) idx[which(refs == ref)[1]]
           else if (any(is.na(refs))) idx[which(is.na(refs))[1]]
           else NA_integer_
    if (is.na(pick)) {
      message("Gene panel '", nm, "' ships only for reference(s) ",
              paste(unique(refs), collapse = ", "), " — not offered for ",
              if (is.na(ref)) "an unknown reference" else ref)
      next
    }
    p = tryCatch(load_gene_panel(tsv_files[pick], reference),
      error = function(e) {
        message("Skipping gene panel '", nm, "': ", conditionMessage(e)); NULL })
    if (!is.null(p)) panels[[nm]] = p
  }
  panels
}

# ---- Small JS serialisation helpers --------------------------------------
# The report ships panel data and column positions to its own client-side filter.
# Keeping these here means the R table and the JS that indexes it are generated
# from the same object, instead of the JS re-deriving positions from the rendered
# header (which breaks under DT's filter row and Scroller's cloned thead).

js_quote = function(x) paste0('"', gsub('"', '\\\\"', as.character(x)), '"')

# {"col":0,"other":1} — column name to zero-based index, for a rownames=FALSE DT.
js_col_index_map = function(nms) {
  if (length(nms) == 0) return("{}")
  paste0("{", paste0(js_quote(nms), ":", seq_along(nms) - 1L, collapse = ","), "}")
}

# One JS number literal per element. Element-wise rather than vectorised because
# format() on a vector pads every element to a common format — c(1, 2.5) becomes
# "1.0","2.5" — and scientific notation on a base-pair coordinate would read back as a
# different number than the one it names.
js_num = function(x) {
  vapply(x, function(v) {
    if (is.na(v)) "null" else format(v, scientific = FALSE, trim = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

.js_cell = function(v) if (is.numeric(v)) js_num(v) else ifelse(is.na(v), "null", js_quote(v))

# A JS array literal from an atomic vector. Numbers are emitted bare, everything else
# quoted; NA becomes null so the client can test for it. Written by hand rather than with
# jsonlite because that would be a new dependency in both recipe/meta.yaml and the
# pipeline's environment.yml — see the "R package requirements" note in CLAUDE.md.
js_vec = function(x) {
  if (length(x) == 0) return("[]")
  paste0("[", paste(.js_cell(x), collapse = ","), "]")
}

# A JS array of arrays, one inner array per row of `dt`, columns in `cols` order.
# Row-major and positional: far smaller than an array of objects, which matters when the
# payload is a few thousand cytobands inlined into a self-contained HTML file.
js_rows = function(dt, cols) {
  if (is.null(dt) || nrow(dt) == 0) return("[]")
  cells = lapply(cols, function(cl) .js_cell(dt[[cl]]))
  rows = do.call(paste, c(cells, sep = ","))
  paste0("[[", paste(rows, collapse = "],["), "]]")
}

# Format a number for human-readable display
fmt_bp = function(x) {
  x = as.numeric(x)
  ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), " Mb"),
    ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 1), " kb"),
      paste0(x, " bp")))
}

# Embed a local PNG file as a self-contained base64 img tag
embed_png = function(path, max_width = "900px") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  b64 = base64enc::base64encode(path)
  htmltools::tags$img(
    src   = paste0("data:image/png;base64,", b64),
    style = paste0("max-width:", max_width, "; display:block; margin:auto;")
  )
}

# Build a data: URI for a Wakhan Plotly HTML file, with a small responsive-resize
# script injected before </body> so the plot fills the iframe's width instead of
# rendering at Plotly's fixed native layout.width (which causes horizontal scroll
# inside the iframe). Runs on the iframe's own `load` event so it fires after
# Plotly.newPlot() has already drawn the figure.
wakhan_plot_datauri = function(path) {
  html = paste(readLines(path, warn = FALSE), collapse = "\n")
  # Wakhan's Plotly divs carry an inline fixed width/height (e.g. style=\"width:1380px\")
  # set by Plotly at export time, in addition to a fixed layout.width. autosize/relayout
  # alone resizes against that fixed div, so the div's own inline size must be cleared
  # to 100% first, then relayout({autosize:true}) + Plots.resize() recomputes against
  # the now-flexible container (i.e. the iframe).
  resize_script = "
<script>
window.addEventListener('load', function () {
  function resizeAll() {
    if (!window.Plotly) return;
    document.querySelectorAll('.plotly-graph-div').forEach(function (gd) {
      gd.style.width = '100%';
      gd.style.height = '100%';
      Plotly.relayout(gd, {autosize: true});
      Plotly.Plots.resize(gd);
    });
  }
  resizeAll();
  window.addEventListener('resize', resizeAll);
});
</script>
"
  if (grepl("</body>", html, fixed = TRUE)) {
    html = sub("</body>", paste0(resize_script, "</body>"), html, fixed = TRUE)
  } else {
    html = paste0(html, resize_script)
  }
  paste0("data:text/html;base64,", base64enc::base64encode(charToRaw(html)))
}

# Embed a self-contained HTML file (e.g. a standalone Plotly plot) as an inline iframe.
embed_html_iframe = function(path, height = "780px") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  htmltools::tags$iframe(
    src   = wakhan_plot_datauri(path),
    style = paste0("width:100%; height:", height, "; border:none;")
  )
}

# Render Wakhan's ranked copy-number plots as a self-contained tab widget (not a
# Quarto .panel-tabset): Quarto's panel-tabset relies on Pandoc parsing `####`
# ATX headings out of a results='asis' stream, which breaks when raw iframe HTML
# for one rank is emitted immediately before the next rank's heading (Pandoc
# absorbs the heading into the preceding raw-HTML block, so only the first tab
# ever registers). This widget also defers loading: only the first pane's
# iframe gets a real `src`; the rest carry `data-src` and are populated on
# first click, so hidden ranks' plotly.js payloads aren't parsed at page load.
render_wakhan_cn_tabs = function(plots) {
  if (length(plots) == 0) return(NULL)

  ids = paste0("wakhan-cn-pane-", seq_along(plots))

  buttons = lapply(seq_along(plots), function(i) {
    p = plots[[i]]
    htmltools::tags$button(
      class = if (i == 1) "wakhan-cn-tab active" else "wakhan-cn-tab",
      `data-target` = ids[i],
      paste0("Rank ", p$rank, " — purity ", p$purity, ", ploidy ", p$ploidy)
    )
  })

  panes = lapply(seq_along(plots), function(i) {
    p = plots[[i]]
    uri = wakhan_plot_datauri(p$plot)
    iframe = if (i == 1) {
      htmltools::tags$iframe(src = uri, style = "width:100%; height:780px; border:none;")
    } else {
      htmltools::tags$iframe(`data-src` = uri, style = "width:100%; height:780px; border:none;")
    }
    htmltools::tags$div(
      class = if (i == 1) "wakhan-cn-pane active" else "wakhan-cn-pane",
      id    = ids[i],
      iframe
    )
  })

  htmltools::tagList(
    htmltools::tags$style("
      .wakhan-cn-tabs__nav { display:flex; flex-wrap:wrap; gap:6px; margin-bottom:10px; }
      .wakhan-cn-tab {
        border:1px solid var(--color-border, #ccc); background:var(--color-bg, #fff);
        border-radius:5px; padding:5px 10px; font-size:0.85rem; cursor:pointer;
      }
      .wakhan-cn-tab.active { background:var(--color-primary, #333); color:#fff; }
      .wakhan-cn-pane { display:none; }
      .wakhan-cn-pane.active { display:block; }
    "),
    htmltools::tags$div(class = "wakhan-cn-tabs__nav", buttons),
    htmltools::tags$div(class = "wakhan-cn-tabs__panes", panes),
    htmltools::tags$script(htmltools::HTML("
      document.querySelectorAll('.wakhan-cn-tabs__nav').forEach(function (nav) {
        nav.querySelectorAll('.wakhan-cn-tab').forEach(function (btn) {
          btn.addEventListener('click', function () {
            const container = nav.nextElementSibling;
            nav.querySelectorAll('.wakhan-cn-tab').forEach(function (b) { b.classList.remove('active'); });
            btn.classList.add('active');
            container.querySelectorAll('.wakhan-cn-pane').forEach(function (p) { p.classList.remove('active'); });
            const pane = document.getElementById(btn.dataset.target);
            pane.classList.add('active');
            const iframe = pane.querySelector('iframe[data-src]');
            if (iframe) {
              iframe.src = iframe.dataset.src;
              iframe.removeAttribute('data-src');
            }
          });
        });
      });
    "))
  )
}

# Compute coding TMB from a variant_table produced by build_variant_table().
# consequence column may be comma-joined (e.g. "frameshift_variant,splice_region_variant").
# denominator_mb: coding Mb used as divisor (default 30 Mb — canonical clinical denominator).
compute_tmb = function(variant_table, denominator_mb = 30) {
  nonsyn_terms = c(
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant"
  )
  if (is.null(variant_table) || nrow(variant_table) == 0) {
    return(list(n_nonsyn = NA_integer_, tmb = NA_real_, denominator_mb = denominator_mb))
  }
  is_nonsyn = vapply(variant_table$consequence, function(csq) {
    if (is.na(csq) || csq == "") return(FALSE)
    any(trimws(unlist(strsplit(csq, ","))) %in% nonsyn_terms)
  }, logical(1))
  n_nonsyn = sum(is_nonsyn, na.rm = TRUE)
  list(
    n_nonsyn       = n_nonsyn,
    tmb            = round(n_nonsyn / denominator_mb, 2),
    denominator_mb = denominator_mb
  )
}
