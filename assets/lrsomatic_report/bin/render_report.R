#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(quarto)
  library(yaml)
})

# Locate the repository root relative to this script. normalizePath() must be
# applied to the script *file* path (resolving a bin/ symlink to its real
# target) before taking dirname() -- doing it the other way around resolves
# the symlink's containing directory instead, which is a no-op when that
# directory isn't itself a symlink (e.g. $PREFIX/bin from the bioconda recipe).
script_file = normalizePath(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
repo_dir    = normalizePath(file.path(dirname(script_file), ".."))

# Source helpers (needed for detect_reference and locate_outputs below)
source(file.path(repo_dir, "R/utils.R"))
source(file.path(repo_dir, "R/references.R"))
source(file.path(repo_dir, "R/locate_outputs.R"))

# ---- CLI argument parsing -----------------------------------------------
option_list = list(
  make_option("--sample-dir",  type = "character", default = NULL,
              help = "Path to the sample output directory (required)"),
  make_option("--sample-id",   type = "character", default = NULL,
              help = "Sample identifier, e.g. SAMPLE_ID (required)"),
  make_option("--reference",   type = "character", default = "auto",
              help = "Reference genome: t2t | hg38 | auto (default: auto)"),
  make_option("--sex",         type = "character", default = NULL,
              help = "Biological sex: male | female | XY | XX (required)"),
  make_option("--gene-panel",  type = "character", default = "none",
              help = paste("Gene panel applied on load: none | builtin name (lymphoid) | path to TSV",
                           "(default: none, i.e. unfiltered). Repeatable — pass it several times to",
                           "open with several panels applied at once; a variant or SV is kept if it",
                           "hits any of them.")),
  make_option("--output",      type = "character", default = NULL,
              help = "Output HTML path (default: <sample-id>_report.html in current dir)"),
  make_option("--title",       type = "character", default = NULL,
              help = "Report title (default: 'LRSomatic Report – <sample-id>')")
)

# --gene-panel is repeatable, which optparse cannot express (it has no action="append"
# and would keep only the last value). Strip every occurrence from argv first and hand
# parse_args() the remainder; --help and --version are untouched by the pre-scan.
argv        = commandArgs(trailingOnly = TRUE)
gene_panel_args = extract_repeated_option(argv, "--gene-panel")
opt         = parse_args(OptionParser(option_list = option_list), args = gene_panel_args$rest)
gene_panels = if (length(gene_panel_args$values) == 0) "none" else gene_panel_args$values

# ---- Validate required arguments ----------------------------------------
abort = function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }

if (is.null(opt[["sample-dir"]])) abort("--sample-dir is required")
if (is.null(opt[["sex"]]))        abort("--sex is required")

sample_dir  = normalizePath(opt[["sample-dir"]], mustWork = TRUE)
sample_id   = if (!is.null(opt[["sample-id"]])) opt[["sample-id"]] else basename(sample_dir)
sex         = tolower(trimws(opt[["sex"]]))
sex         = switch(sex, xy = "male", xx = "female", sex)  # normalise XY/XX

output     = if (!is.null(opt[["output"]])) opt[["output"]] else
             file.path(getwd(), paste0(sample_id, "_report.html"))
title      = if (!is.null(opt[["title"]])) opt[["title"]] else
             paste0("LRSomatic Report – ", sample_id)

# ---- Locate per-tool outputs ---------------------------------------------
message("Locating outputs in: ", sample_dir)
outputs = locate_outputs(sample_dir, sample_id)
message("Run mode: ", outputs$mode)
message("VEP somatic: ", ifelse(is.null(outputs$vep_somatic), "NOT FOUND", outputs$vep_somatic))
message("Somatic VAF VCF: ", ifelse(is.null(outputs$somatic_vaf_vcf), "NOT FOUND",
                                    paste(outputs$somatic_vaf_vcf, collapse = ", ")))
message("ASCAT segments: ", ifelse(is.null(outputs$ascat_segments), "NOT FOUND", outputs$ascat_segments))

# ---- Auto-detect reference -----------------------------------------------
# Resolved before the gene panels below, which are reference-specific: a panel's
# coordinates are only valid for the genome they were built on.
reference = opt[["reference"]]
if (reference == "auto") {
  # Reuse already-resolved paths rather than a fixed vep/somatic/* glob
  vep_file = outputs$vep_somatic
  sv_file  = if (is.null(vep_file)) {
    hits = list.files(sample_dir, pattern = "severus_somatic\\.vcf\\.gz$", recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) hits[1] else NA_character_
  } else NA_character_
  probe    = if (!is.null(vep_file)) vep_file else if (!is.na(sv_file)) sv_file else NA_character_
  reference = if (!is.na(probe)) detect_reference(probe) else "t2t"
  message("Auto-detected reference: ", reference)
}
reference = tolower(reference)

# ---- Load all available gene panels ----------------------------------------
# The rendered report always ships every builtin panel so the reader can switch
# panels client-side; --gene-panel only decides which ones are checked on load.
# "__all__" is the sentinel for "nothing selected"; it survives here and in the
# params default of templates/per_sample.qmd, but no longer reaches the browser —
# in the report, no panel checked *is* the unfiltered state.
# Builtins that ship per reference ("lymphoid.hg38.tsv") resolve to one entry for
# the reference detected above; a coordinate panel declaring a different one is a
# hard error rather than a filter matching the wrong genome.
assets_dir = file.path(repo_dir, "assets")
all_panels = load_all_gene_panels(assets_dir, reference)

# "none" means unfiltered, so combining it with a real panel is contradictory rather
# than a case where one of the two quietly wins.
if (any(vapply(gene_panels, is_no_gene_panel, logical(1))) && length(gene_panels) > 1) {
  abort("--gene-panel none cannot be combined with other panels; drop the 'none'.")
}
# Deduplicate by resolved path so the same TSV passed twice registers once rather
# than a second time under a "-custom" key.
gene_panels = unique(vapply(gene_panels, function(g)
  if (file.exists(g)) normalizePath(g) else g, character(1)))

default_panels = character(0)
for (gp in gene_panels) {
  if (is_no_gene_panel(gp)) next
  if (!is.null(builtin_panel_path(assets_dir, gp, reference))) {
    # Load it here too, so a builtin that fails to resolve against this reference aborts
    # now rather than part-way through the Quarto render.
    invisible(tryCatch(resolve_gene_panel(gp, assets_dir, reference),
                       error = function(e) abort(conditionMessage(e))))
    default_panels = c(default_panels, gp)
  } else if (file.exists(gp)) {
    # A user-supplied TSV: register it alongside the builtins so it can be
    # selected on load (and switched away from and back to) in the report.
    nm = unique_panel_name(tools::file_path_sans_ext(basename(gp)), names(all_panels))
    all_panels[[nm]] = tryCatch(load_gene_panel(gp, reference),
                                error = function(e) abort(conditionMessage(e)))
    default_panels = c(default_panels, nm)
  } else {
    abort(paste0("--gene-panel not found: tried builtin '", gp,
                 "' and as a file path. Use 'none' for no filtering."))
  }
}
default_panels = unique(default_panels)
# Kept as a sentinel rather than character(0): an empty vector round-trips through
# Quarto's YAML execute_params as NULL, not as an empty character vector.
if (length(default_panels) == 0) default_panels = "__all__"
message("Gene panels selected on load: ", paste(default_panels, collapse = ", "))

# ---- Render the Quarto template -----------------------------------------
# Copy templates/ and assets/ into a writable working directory: repo_dir's
# own templates/ may be read-only (e.g. inside a container), and Quarto
# writes intermediate files next to the .qmd during render.
work = file.path(getwd(), "._render")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
invisible(file.copy(file.path(repo_dir, "templates"), work, recursive = TRUE))
invisible(file.copy(file.path(repo_dir, "assets"), work, recursive = TRUE))
template = file.path(work, "templates", "per_sample.qmd")
if (!file.exists(template)) abort("Quarto template not found: ", template)

message("Rendering report to: ", output)
quarto::quarto_render(
  input          = template,
  output_file    = basename(output),
  output_format  = "html",
  execute_params = list(
    sample_id     = sample_id,
    sample_dir    = sample_dir,
    reference     = reference,
    sex           = sex,
    default_panels = default_panels,
    all_panels    = all_panels,
    title         = title,
    repo_dir      = repo_dir,
    outputs       = outputs
  ),
  quiet = FALSE
)

# Move output if Quarto wrote it next to the template
rendered = file.path(dirname(template), basename(output))
if (file.exists(rendered)) {
  dest = normalizePath(output, mustWork = FALSE)
  src  = normalizePath(rendered, mustWork = FALSE)
  if (src != dest) {
    ok = file.copy(rendered, output, overwrite = TRUE)
    if (ok) file.remove(rendered)
  }
}
unlink(work, recursive = TRUE)

if (file.exists(output)) {
  message("Report written to: ", output)
} else {
  abort("Rendering completed but output file not found at: ", output)
}
