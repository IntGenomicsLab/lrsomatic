# Structural variants section. Reference implementation of the section-module
# contract (see R/sections.R and CLAUDE.md). Keyed by caller so a second SV
# caller can be added later without touching the plumbing below.

register_section(list(
  id = "sv",
  title = "Structural variants",

  locate = function(sample_dir, sample_id) {
    d = sample_dir

    find1 = function(pattern) {
      hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) hits[1] else NULL
    }

    severus_vcf      = find1("^severus_somatic\\.vcf\\.gz$")
    severus_gene_tsv = find1("^SV_filtered_with_gene_annotations\\.tsv$")
    # VEP SV VCF (CSQ-annotated) is the more commonly produced annotation source; the
    # gene-annotated TSV above is a fallback for pipelines that produce it instead.
    severus_vep_vcf  = find1("_SV_VEP\\.vcf\\.gz$")

    list(callers = list(
      severus = list(vcf = severus_vcf, gene_tsv = severus_gene_tsv, vep_vcf = severus_vep_vcf)
    ))
  },

  parse = function(inputs, section_data) {
    tabs = list()
    circ = list(nontrans = data.table(), translocations = data.table())
    any_annotation = FALSE

    for (nm in names(inputs$callers)) {
      caller_inputs = inputs$callers[[nm]]
      v = parse_severus_vcf(caller_inputs$vcf)

      # VEP SV VCF is the primary annotation source; the gene-annotated TSV (not produced
      # by most pipelines) is a fallback for samples that have it instead.
      if (!is.null(caller_inputs$vep_vcf)) {
        t = build_sv_table_from_vep(caller_inputs$vcf, caller_inputs$vep_vcf)
      } else {
        g = parse_severus_gene_tsv(caller_inputs$gene_tsv)
        t = build_sv_table(g, gene_panel = NULL)
      }
      if (!is.null(t) && nrow(t) > 0) {
        any_annotation = TRUE
        t[, caller := nm]
        tabs[[nm]] = t
      }
      # Circos tracks are drawn from raw breakpoints, not the gene table;
      # with a single caller today, last-write-wins is a no-op.
      circ$nontrans = v$nontrans
      circ$translocations = v$translocations
    }

    tbl = if (length(tabs) > 0) rbindlist(tabs, fill = TRUE) else data.table()

    list(
      table            = tbl,
      circos           = circ,
      annotation_found = any_annotation
    )
  }
))
