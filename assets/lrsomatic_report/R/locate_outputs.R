# Discover the per-tool output files for a sample run. Discovery is recursive
# under sample_dir: the pipeline may dump inputs flat rather than in a fixed
# directory tree, so files are matched by their distinctive filename suffix.
# Returns a named list; any missing optional file is NULL.

locate_outputs = function(sample_dir, sample_id) {
  d = sample_dir  # shorthand

  # First recursive hit under `root` matching a filename pattern
  find1 = function(pattern, root = d) {
    hits = list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) hits[1] else NULL
  }

  # Same, but excluding anything under a normal/ subtree (tumor-side QC)
  find1_tumor = function(pattern) {
    hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
    hits = hits[!grepl("/normal/", hits)]
    if (length(hits) > 0) hits[1] else NULL
  }

  # The mirror of find1_tumor: normal-side files, wherever the pipeline puts them.
  # This has moved (a top-level normal/ historically, qc/normal/ today), so match on
  # the path component rather than rooting the search at a fixed directory.
  find1_normal = function(pattern) {
    hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
    hits = hits[grepl("/normal/", hits)]
    if (length(hits) > 0) hits[1] else NULL
  }

  # --- small variants -------------------------------------------------------
  vep_somatic = find1("_SOMATIC_VEP\\.vcf\\.gz$")

  # VAF, depth and phasing come from the VCF that VEP annotated, not from a separate
  # caller VCF. Preferred is the phased somatic VCF, which is what the pipeline feeds to
  # VEP and which additionally carries FORMAT/PS. Runs predating variants/phased/ fall back
  # to the raw ClairS(-TO) output, matched on the containing directory because the
  # basename ("somatic.vcf.gz") is shared across callers. The clairs*/ fallback may return
  # several paths — parse_caller_vcf() stacks them.
  somatic_vaf_vcf = {
    phased = file.path(d, "variants", "phased", "somatic_smallvariants.vcf.gz")
    if (file.exists(phased)) phased else {
      vcfs = list.files(d, pattern = "\\.vcf\\.gz$", recursive = TRUE, full.names = TRUE)
      named = vcfs[grepl("/clairs(to)?/somatic\\.vcf\\.gz$", vcfs)]
      any_c = vcfs[grepl("/clairs(to)?/", vcfs) & !grepl("/germline\\.vcf\\.gz$", vcfs)]
      if (length(named) > 0) named[1] else if (length(any_c) > 0) any_c else NULL
    }
  }

  # --- structural variants ---------------------------------------------------
  # Severus paths are now located by R/sections/sv.R (section-module contract).

  # --- ASCAT ----------------------------------------------------------------
  ascat_segments_raw = find1("\\.segments_raw\\.txt$")
  ascat_purityploidy = find1("\\.purityploidy\\.txt$")
  ascat_plots = list(
    profile    = find1("\\.tumour\\.ASCATprofile\\.png$"),
    rawprofile = find1("\\.tumour\\.rawprofile\\.png$"),
    sunrise    = find1("\\.tumour\\.sunrise\\.png$"),
    aspcf      = find1("\\.tumour\\.ASPCF\\.png$"),
    before_gc  = find1("\\.before_correction\\..*\\.tumour\\.tumour\\.png$"),
    after_gc   = find1("\\.after_correction_gc.*\\.tumour\\.tumour\\.png$"),
    tumour_sep = find1("^tumorSep.*\\.tumour\\.png$")
  )

  # --- QC (tumor side) --------------------------------------------------------
  mosdepth_summary = find1_tumor("\\.mosdepth\\.summary\\.txt$")
  mosdepth_dist    = find1_tumor("\\.mosdepth\\.global\\.dist\\.txt$")
  cramino_aln      = find1_tumor("_cramino\\.txt$")
  flagstat         = find1_tumor("\\.flagstat$")
  samtools_stats   = find1_tumor("\\.stats$")

  # --- Normal-side QC (matched mode only) -----------------------------------
  normal_mosdepth_summary = find1_normal("\\.mosdepth\\.summary\\.txt$")
  normal_mosdepth_dist    = find1_normal("\\.mosdepth\\.global\\.dist\\.txt$")
  normal_cramino          = find1_normal("_cramino\\.txt$")
  normal_flagstat         = find1_normal("\\.flagstat$")
  normal_samtools_stats   = find1_normal("\\.stats$")

  # Driven by what was actually found rather than by directory layout or run mode:
  # the QC section renders a tumour/normal comparison only if there is normal data.
  has_normal = !is.null(normal_mosdepth_summary) || !is.null(normal_cramino)

  # Run mode is derived from the same evidence rather than declared by the caller; its
  # only consumer is the hero badge in templates/sections/_header.qmd.
  mode = if (has_normal) "matched" else "tumour-only"

  # --- Wakhan (optional) -----------------------------------------------------
  wakhan_dir = file.path(d, "wakhan")
  has_wakhan = dir.exists(wakhan_dir)
  wakhan_solutions = if (has_wakhan) {
    f = file.path(wakhan_dir, "solutions_ranks.tsv")
    if (file.exists(f)) f else NULL
  } else NULL
  wakhan_heatmap = if (has_wakhan) {
    hits = Sys.glob(file.path(wakhan_dir, "*heatmap_ploidy_purity.html"))
    if (length(hits) > 0) hits[1] else NULL
  } else NULL

  list(
    mode             = mode,
    vep_somatic      = vep_somatic,
    somatic_vaf_vcf  = somatic_vaf_vcf,
    ascat_segments   = ascat_segments_raw,
    ascat_purityploidy = ascat_purityploidy,
    mosdepth_summary = mosdepth_summary,
    mosdepth_dist    = mosdepth_dist,
    cramino          = cramino_aln,
    flagstat         = flagstat,
    samtools_stats   = samtools_stats,
    has_normal       = has_normal,
    ascat_plots      = ascat_plots,
    normal_mosdepth_summary = normal_mosdepth_summary,
    normal_mosdepth_dist    = normal_mosdepth_dist,
    normal_cramino          = normal_cramino,
    normal_flagstat         = normal_flagstat,
    normal_samtools_stats   = normal_samtools_stats,
    has_wakhan       = has_wakhan,
    wakhan_dir       = wakhan_dir,
    wakhan_solutions = wakhan_solutions,
    wakhan_heatmap   = wakhan_heatmap
  )
}
