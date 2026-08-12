suppressPackageStartupMessages({
  library(data.table)
})

ensure_chr_prefix = function(x) {
  ifelse(startsWith(x, "chr"), x, paste0("chr", x))
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

# Load a gene panel TSV or plain text file; returns character vector of gene symbols
load_gene_panel = function(path) {
  if (!file.exists(path)) stop("Gene panel file not found: ", path)
  dt = tryCatch(
    fread(path, header = TRUE, sep = "\t", fill = TRUE),
    error = function(e) fread(path, header = FALSE, sep = "\t", fill = TRUE)
  )
  gene_col = if ("gene" %in% tolower(names(dt))) names(dt)[tolower(names(dt)) == "gene"][1] else names(dt)[1]
  unique(dt[[gene_col]])
}

# Filter a data frame to rows where the gene column matches the panel.
# panel_genes = NULL means "no panel" and returns dt untouched; an empty
# character vector is a genuinely empty panel and filters everything out.
filter_by_gene_panel = function(dt, panel_genes, gene_col = "gene") {
  if (is.null(panel_genes)) return(dt)
  dt[dt[[gene_col]] %in% panel_genes, ]
}

# Is a --gene-panel argument the "no filtering" sentinel?
is_no_gene_panel = function(panel_arg) {
  is.null(panel_arg) || length(panel_arg) != 1 || is.na(panel_arg) ||
    identical(tolower(trimws(panel_arg)), "none")
}

# Resolve a --gene-panel arg: the "none" sentinel (no filtering, returns NULL),
# a builtin name ("lymphoid"), or a path to a TSV. A value that is neither is an
# error rather than a silent fall-through to unfiltered output.
resolve_gene_panel = function(panel_arg, assets_dir) {
  if (is_no_gene_panel(panel_arg)) return(NULL)
  builtin_path = file.path(assets_dir, "gene_lists", paste0(panel_arg, ".tsv"))
  if (file.exists(builtin_path)) return(load_gene_panel(builtin_path))
  if (file.exists(panel_arg)) return(load_gene_panel(panel_arg))
  stop("Gene panel not found (tried builtin '", panel_arg, "' and as file path)")
}

# Load all gene panels from assets/gene_lists/*.tsv
# Returns a named list (name = panel name, value = character vector of gene symbols)
load_all_gene_panels = function(assets_dir) {
  tsv_files = Sys.glob(file.path(assets_dir, "gene_lists", "*.tsv"))
  if (length(tsv_files) == 0) return(list())
  panels = lapply(tsv_files, load_gene_panel)
  names(panels) = tools::file_path_sans_ext(basename(tsv_files))
  panels
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
