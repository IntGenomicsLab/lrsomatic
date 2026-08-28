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
              help = "Gene panel applied on load: none | builtin name (lymphoid) | path to TSV (default: none, i.e. unfiltered)"),
  make_option("--output",      type = "character", default = NULL,
              help = "Output HTML path (default: <sample-id>_report.html in current dir)"),
  make_option("--title",       type = "character", default = NULL,
              help = "Report title (default: 'LRSomatic Report – <sample-id>')")
)

opt = parse_args(OptionParser(option_list = option_list))

# ---- Validate required arguments ----------------------------------------
abort = function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }

if (is.null(opt[["sample-dir"]])) abort("--sample-dir is required")
if (is.null(opt[["sex"]]))        abort("--sex is required")

sample_dir  = normalizePath(opt[["sample-dir"]], mustWork = TRUE)
sample_id   = if (!is.null(opt[["sample-id"]])) opt[["sample-id"]] else basename(sample_dir)
sex         = tolower(trimws(opt[["sex"]]))
sex         = switch(sex, xy = "male", xx = "female", sex)  # normalise XY/XX

gene_panel = opt[["gene-panel"]]
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
# panels client-side; --gene-panel only decides which one is selected on load.
# "__all__" is the sentinel the report's JS uses for "no filter" — it must stay
# in sync with templates/sections/_gene_filter.qmd and the search hook in
# templates/per_sample.qmd.
# Builtins that ship per reference ("lymphoid.hg38.tsv") resolve to one entry for
# the reference detected above; a coordinate panel declaring a different one is a
# hard error rather than a filter matching the wrong genome.
all_panels = load_all_gene_panels(file.path(repo_dir, "assets"), reference)
default_panel = if (is_no_gene_panel(gene_panel)) {
  gene_panel = "none"
  "__all__"
} else if (!is.null(builtin_panel_path(file.path(repo_dir, "assets"), gene_panel, reference))) {
  # Load it here too, so a builtin that fails to resolve against this reference aborts
  # now rather than part-way through the Quarto render.
  invisible(tryCatch(resolve_gene_panel(gene_panel, file.path(repo_dir, "assets"), reference),
                     error = function(e) abort(conditionMessage(e))))
  gene_panel
} else if (file.exists(gene_panel)) {
  # A user-supplied TSV: register it alongside the builtins so it can be
  # selected on load (and switched away from and back to) in the report.
  nm = tools::file_path_sans_ext(basename(gene_panel))
  if (nm %in% names(all_panels)) nm = paste0(nm, "-custom")
  all_panels[[nm]] = tryCatch(load_gene_panel(gene_panel, reference),
                              error = function(e) abort(conditionMessage(e)))
  # Absolute, because the template resolves it again from Quarto's own working
  # directory (the copied template dir), not from where this script was invoked.
  gene_panel = normalizePath(gene_panel)
  nm
} else {
  abort(paste0("--gene-panel not found: tried builtin '", gene_panel,
               "' and as a file path. Use 'none' for no filtering."))
}

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
    gene_panel    = gene_panel,
    default_panel = default_panel,
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
