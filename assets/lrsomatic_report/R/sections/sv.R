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
    annotation_path = NULL

    for (nm in names(inputs$callers)) {
      caller_inputs = inputs$callers[[nm]]

      # One parse of the caller VCF feeds both the table and the circos tracks, so the
      # SV count, the panel filter and the drawn links cannot disagree about what a
      # rearrangement is.
      records = parse_severus_somatic_records(caller_inputs$vcf)

      if (nrow(records) > 0) {
        # The VEP SV VCF only supplies the per-breakend display symbols; a missing one
        # leaves those columns empty rather than losing the SVs, because panel matching
        # is done on coordinates (see sv_panel_hits()).
        t = build_sv_table_from_vep(records, caller_inputs$vep_vcf)
        if (!is.null(caller_inputs$vep_vcf)) annotation_path = caller_inputs$vep_vcf
      } else {
        # No caller VCF: the gene-annotated TSV is the fallback for pipelines that
        # publish it instead.
        t = build_sv_table(parse_severus_gene_tsv(caller_inputs$gene_tsv))
        if (nrow(t) > 0) annotation_path = caller_inputs$gene_tsv
      }

      if (!is.null(t) && nrow(t) > 0) {
        t[, caller := nm]
        tabs[[nm]] = t
      }
      # Circos tracks are drawn from the same mate-collapsed records as the table;
      # with a single caller today, last-write-wins is a no-op.
      circ = severus_circos_tracks(records)
    }

    tbl = if (length(tabs) > 0) rbindlist(tabs, fill = TRUE) else data.table()

    list(
      table            = tbl,
      circos           = circ,
      annotation_path  = annotation_path,
      annotation_found = !is.null(annotation_path)
    )
  }
))
