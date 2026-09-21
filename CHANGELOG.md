# IntGenomicsLab/lrsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### `Added`

- [#190](https://github.com/IntGenomicsLab/lrsomatic/pull/190) - Added mutational signature analysis: SigProfilerMatrixGenerator matrices (SBS/DBS/ID) and per-sample COSMIC signature fitting with SigProfilerAssignment on the phased somatic VCF, with CHM13-T2T support via a purpose-built image (`ghcr.io/ljwharbers/sigprofiler`) and CHM13-renormalised COSMIC signatures; new `--sigprofiler_*` / `--download_sigprofiler_genome` / `--skip_signatures` parameters (@ljwharbers).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - Added `modkit_args` (default `--cpg --modified-bases 5mC`) to control the arguments passed to `modkit pileup`, and `modkit_phased` to run the pileup with `--phased` on the Longphase-haplotagged BAMs, producing `_hp1`, `_hp2` and `_combined` bedMethyl files per sample (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added `LRSOMATICREPORT` as the final pipeline step: a self-contained per-sample HTML report covering small variants, structural variants, copy number and QC. Skip it with `--skip_report`; choose the gene panel selected on load with `--report_gene_panel` (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Vendored the [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) v1.3.2 tool source at `assets/lrsomatic_report`, so `nextflow run IntGenomicsLab/lrsomatic` ships it without a submodule checkout (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added a `solution_dirs` output to the WAKHAN module so its per-solution copy-number plots can be staged downstream (@ljwharbers).
- [#193](https://github.com/IntGenomicsLab/lrsomatic/pull/193) - Added VEP plugins for AlphaMissense, SIFT/PolyPhen, ClinVar and REVEL on GRCh38, and AlphaMissense plus SIFT/PolyPhen on CHM13 via protein-space lookup. Enabled by default with `--genome GRCh38` or `--genome CHM13`: the pipeline fetches each resource from its original source and prepares the releases VEP cannot read as published, so a first GRCh38 run downloads around 1.4 GB. Any `--vep_*` path overrides its default, and `--skip_vep_plugins` turns the set off. CADD and EVE are wired up but left opt-in, behind `--vep_cadd_snv`/`--vep_cadd_indel` and `--vep_eve`, since they are 81 GB and 9.6 GB respectively. Every default points at the resource's original source, except the AlphaMissense GRCh38 tabix index and the CHM13 protein-space table, which the lab hosts and which are CC BY 4.0 with attribution in `CITATIONS.md`. Prepared REVEL and EVE files are published to `<outdir>/vep_plugins/` so later runs can reuse them (@AmberVerhasselt).

### `Changed`

- [#XXX](https://github.com/IntGenomicsLab/lrsomatic/pull/XXX) - The somatic arm is now recovered from the phased germline+somatic VCF by provenance rather than by position. `PHASING_HAPLOTYPING:BCFTOOLS_VIEW` previously used the somatic VCF as a `-T` targets file, which matches on `CHROM`/`POS` only, so every germline record sitting at a somatic coordinate was retained; on B2037620 it removed none of the 13,708,100 germline records, and `variants/phased/somatic_smallvariants.vcf.gz` ended up holding 4,129,270 germline variants that the report could not distinguish from somatic ones. A new `VCFTAG` module (aliased `TAG_SOMATIC`/`TAG_GERMLINE`) now stamps `INFO/SOMATIC` and `INFO/GERMLINE` on the two arms immediately before they are merged for phasing — the only point at which origin is unambiguous for every caller, since `GERMLINE_CONSENSUS` can emit records that never passed through `VCFSPLIT` — and the somatic arm is selected with `-i 'INFO/SOMATIC=1'`. Longphase preserves custom INFO keys, verified on v2.0.1. `VCFSPLIT` additionally records each record's original `FILTER` in `INFO/ORIG_FILTER` on both splits (multi-valued `FILTER` is joined with `,`, since `;` separates INFO fields); it still normalises the germline split's `FILTER` to `PASS` so that downstream tools which filter on `PASS` continue to see every record, but that value is no longer destroyed. Germline calls remain published in full under `variants/phased/germline_smallvariants.vcf.gz`, `vep/germline/` and `variants/<caller>/` (@robert-a-forsyth).
- [#XXX](https://github.com/IntGenomicsLab/lrsomatic/pull/XXX) - `*_var_combine = 'all'` now produces the union it is documented to produce. Both branches of `SMALL_VARIANT_CONSENSUS` concatenated the shared record from the prioritized caller with the private calls of the *other* caller, so the prioritized caller's own private calls were always discarded: with the default `prioritize_caller_somatic = 'clair'` that silently dropped ClairS-TO's private calls (99 on B2037620), and with `'deepsomatic'` it would have dropped every DeepSomatic-private call (26,560 `PASS` records on the same sample). Both private sets are now kept alongside the shared record, matching `nextflow_schema.json` ("keeps all variants from both callers"); `prioritize_caller_*` selects only whose record represents a shared variant. Invalid `combine_method`/`prioritize_caller` values now raise a clear error instead of leaving the output channel undefined (@robert-a-forsyth).
- [#XXX](https://github.com/IntGenomicsLab/lrsomatic/pull/XXX) - Small variant caller output is now restricted to `PASS` records before it is used downstream, controlled by the new `--smallvar_filter_pass` parameter (`true` by default). DeepVariant and DeepSomatic emit a record for every site they evaluate rather than only the variants they call: on a 30x PacBio tumour sample `<sample>_somatic.vcf.gz` held 13,684,025 records of which 49,957 were `PASS` (9,349,614 `RefCall`, 4,011,128 `GERMLINE`, 273,326 `PON`), and `<sample>_germline.vcf.gz` held 13,684,023 records with 4,991,797 `PASS`. Clair3 and ClairS are far less extreme but still kept their `LowQual` and `NonSomatic` records, so with the default `*_var_combine = 'all'` the union was effectively "every site every caller looked at": the published `variants/phased/somatic_smallvariants.vcf.gz` reached 27,368,030 records for ~4.1 M unique SNVs, giving a coding TMB around 406 mut/Mb, and `LRSOMATICREPORT` could not render it at all — its circos and small-variant chunks each exceeded V8's 512 MB string limit, so Quarto failed with `failed to allocate string; buffer exceeds maximum length`. DeepVariant, DeepSomatic, Clair3 and ClairS now each pass through a `PASS`-only copy before the caller consensus, phasing, VEP and the report see them; ClairS-TO is unchanged because `VCFSPLIT` already restricted it to `PASS`, and signature fitting is unchanged because `SIGNATURES_BCFTOOLS_VIEW` already applied its own `PASS` filter. Reruns therefore give different `somatic_smallvariants.vcf.gz` and `germline_smallvariants.vcf.gz` content, and VEP and the report see far fewer variants. Set `--smallvar_filter_pass false` for the previous behaviour. The per-caller VCFs published under `variants/deepvariant`, `variants/deepsomatic`, `variants/clair3` and `variants/clairs` are unchanged and remain unfiltered (@robert-a-forsyth).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - `MODKIT_PILEUP` now runs `modkit pileup` with `--cpg --modified-bases 5mC` by default; previously it ran with no arguments, and that unfiltered pileup (every modification code at every sequence context) produced 30-45 GB per sample. Reruns therefore give different bedMethyl content: only 5mC rows at CpG sites, so 5hmC and 6mA calls and non-CpG positions are no longer reported. Set `--modkit_args=''` (the `=` form; `--modkit_args ''` reaches the pipeline as `true` and is rejected by parameter validation on Nextflow 25) or an empty `modkit_args` entry in a params file to get the previous unfiltered output (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - Re-synced the vendored [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) to v1.3.0, which adds tickbox dropdown filters on the categorical columns of both variant tables and turns the report's gene panel selector into checkboxes (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Re-synced the vendored [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) to v1.3.2: facet dropdown counts follow the active filters, opening a facet menu no longer resets the table's horizontal scroll, a flatter clinical theme, and inline code comments trimmed to one line (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - `--report_gene_panel` now takes a comma-separated list, so several panels can be applied at once: a variant or SV is kept if it hits any of them. Panel values are also validated at launch instead of failing inside the report task (@ljwharbers).
- [#184](https://github.com/IntGenomicsLab/lrsomatic/pull/184) - Replaced the CHM13 Severus panel of normals with the merged 1000 Genomes + ASAP panel (@AmberVerhasselt).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Updated to the nf-core template for [nf-core/tools 4.1.0](https://github.com/nf-core/tools/releases/tag/4.1.0). Adds the `check-added-large-files`, `check-merge-conflict` and `block-pipeline-outdir` pre-commit hooks, a `process_low_memory` resource label, and a split-out `pr-comment.yml` workflow. `conf/igenomes.config` was converted to the template's strict-syntax `params.genomes` map literal, and the MultiQC module was bumped to 1.35 (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Raised the minimum Nextflow version to `25.10.4`, following the template, and bumped the nf-test CI matrix to match (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Filled in the `manifest.contributors` `contribution` fields: @ljwharbers and @robert-a-forsyth as author and maintainer, all other contributors as author. Previously empty, which left the RO-Crate metadata without any authors (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Added @AmberVerhasselt to `manifest.contributors` as a contributor and to the README credits, and expanded @ljwharbers' affiliation to match the other core contributors. RO-Crate metadata updated to match (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Set @laulambr's affiliation to the same three-part VIB/KU Leuven string as the other core contributors (@ljwharbers).

### `Fixed`

- [#193](https://github.com/IntGenomicsLab/lrsomatic/pull/193) - Fixed the documented way of enabling EVE. `--vep_eve https://evemodel.org/api/proteins/bulk/download/` was rejected at launch, because the check for "still needs reshaping" keyed on a `.zip` suffix and that endpoint carries no file extension, so the URL was taken for a finished file and an index demanded for it. The check now asks whether the value already is a prepared (bgzipped) file instead (@AmberVerhasselt).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - `MODKIT_PILEUP` now runs a patched modkit 0.6.4 built from [ljwharbers/modkit@pacbio-conflict-fix](https://github.com/ljwharbers/modkit/tree/pacbio-conflict-fix): `ghcr.io/ljwharbers/modkit:0.6.4-pacbiofix-6e0afa2` under Docker, and the native SIF `oras://ghcr.io/ljwharbers/modkit-sif:0.6.4-pacbiofix-6e0afa2` under Singularity/Apptainer (unless `singularity_pull_docker_container` is set, which falls back to the Docker image). It keeps reads whose 5mC and 5hmC probabilities sum above 1, fixes pileup on PacBio-style MM tags, and honours `--phased` and `--modified-bases` in the general pileup workers that modkit uses for PacBio BAMs with 6mA calls (stock modkit wrote empty `_hp1`/`_hp2` files and an `h` row next to every `m` row for those). Stock modkit 0.4.3-0.6.4 silently dropped 32-65 % of reads from recent PacBio HiFi BAMs (Jasmine >= 26.1.3) and returned empty `--cpg` pileups ([nanoporetech/modkit#612](https://github.com/nanoporetech/modkit/issues/612); fixes proposed upstream in [nanoporetech/modkit#720](https://github.com/nanoporetech/modkit/pull/720)). The image is `linux/amd64` only and Conda is not supported for this step: `MODKIT_PILEUP` stops with an error under `-profile conda`/`mamba`, use `--skip_modkit` there. The module should return to the nf-core biocontainer once a modkit release includes the fix (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - NanoPlot pre- and post-alignment statistics now reach MultiQC. `ch_nanoplot_pre_txt` and `ch_nanoplot_post_txt` were declared empty and mixed into the MultiQC inputs, but never assigned from `NANOPLOT_PRE.out.txt` / `NANOPLOT_POST.out.txt`, so the NanoStat section was silently missing from every report (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - `NANOPLOT_PRE` now includes the replicate in its output prefix, so samples with more than one replicate no longer collapse into a single MultiQC sample (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - samtools stats/flagstat/idxstats and mosdepth outputs are now prefixed `<sample>_<type>` instead of `<sample>`. The tumor and normal files of a matched pair shared a basename, so MultiQC logged `Duplicate sample name found! Overwriting` and reported only one of the two; the general statistics table now has one row per BAM (`<sample>_tumor`, `<sample>_normal`), which the post-alignment NanoStat row joins (@ljwharbers).
- [#194](https://github.com/IntGenomicsLab/lrsomatic/pull/194) - A sample whose `clair3_model` differs from the model implied by its BAM header no longer makes Clair3 run twice on every normal BAM sharing that header model. The model download channel was keyed on the header-derived name while downloading the requested model, so `.unique()` kept two entries under one name, `UNTAR` produced two model directories with the same name, and the by-name combine in `PAIRED_SMALLVAR_GERMLINE` ran `CLAIR3` once per copy with the downstream join keeping whichever finished first; the germline VCF, haplotagged BAMs and WhatsHap statistics of the affected sample differed between runs. The same divergence also broke the opposite case, first reported in [#191](https://github.com/IntGenomicsLab/lrsomatic/pull/191): a basecall model absent from the model map made the header-derived name `null`, so `UNTAR` failed with `mkdir: missing operand` even though the explicit `clair3_model` override had downloaded correctly. Entries are now keyed on the downloaded model (@YannVRB, @ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - Raised MultiQC's `log_filesize_limit` to 500 MB. `samtools stats` output for long reads is 70 MB and more, above the 50 MB default, so MultiQC silently skipped it and the samtools stats section covered only small files (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - Stopped snapshotting the md5 of sample4's merged tumour BAM and its index in the `clair_only` nf-test: `samtools merge` gives the colliding `@PG` IDs of the two replicates a random hex suffix, so neither digest is reproducible. The alignment records are, and are now asserted with `bam().getReadsMD5()` instead (@ljwharbers).
- [#182](https://github.com/IntGenomicsLab/lrsomatic/pull/182) - Added `--vcf` to the default `vep_args` so VEP writes VCF output rather than its default tab-delimited format (@AmberVerhasselt).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Corrected the `github` URL for Laurens Lambrechts in `manifest.contributors`, which was a copy of @MariosEft97's, to @laulambr. RO-Crate metadata updated to match (@ljwharbers).

## v1.1.0 - [2026-04-28]

### `Added`

- [#117](https://github.com/IntGenomicsLab/lrsomatic/pull/117) - Added ASCAT PDF plots to output (@robert-a-forsyth).
- [#126](https://github.com/IntGenomicsLab/lrsomatic/pull/126) - Added ASCAT raw segments txt files to output (@AmberVerhasselt).
- [#135](https://github.com/IntGenomicsLab/lrsomatic/pull/135) - Added `skip_m6a` parameter to allow skipping m6A base modification steps (@robert-a-forsyth).
- [#141](https://github.com/IntGenomicsLab/lrsomatic/pull/141) - Added output of phased variants in separate VCF files for improved downstream analysis (@ljwharbers).
- [#143](https://github.com/IntGenomicsLab/lrsomatic/pull/143) - Added Severus `min_support` parameter and `skip_fibernormal` option (@ljwharbers).
- [#145](https://github.com/IntGenomicsLab/lrsomatic/pull/145) - Integrated MultiQC and nanoplot for comprehensive QC reporting with long-read sequencing metrics (@ljwharbers).
- [#147](https://github.com/IntGenomicsLab/lrsomatic/pull/147) - Implemented whatshap_stats module to generate phase block statistics and phasing quality metrics (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added DeepVariant and DeepSomatic modules for germline and somatic variant calling from long-read sequencing data (@robert-a-forsyth).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added GPU support for Clair3, DeepVariant, and fibertools (@robert-a-forsyth).
- [#150](https://github.com/IntGenomicsLab/lrsomatic/pull/150) - Added Claude GitHub Actions workflows for automated code review and PR assistance (@ljwharbers).
- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Integrated modkit module for long-read base modification detection and analysis (@robert-a-forsyth).
- [#165](https://github.com/IntGenomicsLab/lrsomatic/pull/165) - Added bcftools/view and samtools/merge modules; added extended test suites for union, consensus, clair-only, and deep-only caller modes (@robert-a-forsyth).

### `Changed`

- [#123](https://github.com/IntGenomicsLab/lrsomatic/pull/123) - Updated channel structure (@robert-a-forsyth).
- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Bulk module versions update. Fixed some issues with Wakhan (@ljwharbers).
- [#138](https://github.com/IntGenomicsLab/lrsomatic/pull/138) - Perform QC before merging replicates (@robert-a-forsyth).
- [#140](https://github.com/IntGenomicsLab/lrsomatic/pull/140) - Improved documentation with additional pipeline usage examples and configuration guidance (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Refactored variant calling workflow to support both DeepVariant and existing callers with improved configuration handling (@robert-a-forsyth).
- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Updated container versions and dependencies for modkit and related tools (@robert-a-forsyth).
- [#157](https://github.com/IntGenomicsLab/lrsomatic/pull/157) - Added ASAP Panel of Normals citation to CITATIONS.md (@ljwharbers).
- [#160](https://github.com/IntGenomicsLab/lrsomatic/pull/160) - DeepVariant/DeepSomatic optimization and Panel-of-Normals handling improvements; updated docs; removed Claude workflow files (@robert-a-forsyth).
- [#163](https://github.com/IntGenomicsLab/lrsomatic/pull/163) - Added LongPhase supplementary alignment tag to extended args (@AmberVerhasselt).
- [#164](https://github.com/IntGenomicsLab/lrsomatic/pull/164) - Updated nf-core template components (@robert-a-forsyth).
- [#166](https://github.com/IntGenomicsLab/lrsomatic/pull/166) - Updated Wakhan to v0.4.3 using BioContainers distribution (@robert-a-forsyth).

### `Fixed`

- [#116](https://github.com/IntGenomicsLab/lrsomatic/pull/116) - Corrected ASCAT GC and RT bias correction (@AmberVerhasselt).
- [#118](https://github.com/IntGenomicsLab/lrsomatic/pull/118) - Updated nf-core template components to align with latest pipeline standards (@ljwharbers).
- [#130](https://github.com/IntGenomicsLab/lrsomatic/pull/130) - Fixed path pattern for rephased VCF files (@Tim-Yu).
- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Resolved Nextflow strict syntax compliance issues for compatibility with latest Nextflow versions (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Corrected bcftools and vcfsplit operations for accurate variant filtering and merging (@robert-a-forsyth).
- [#165](https://github.com/IntGenomicsLab/lrsomatic/pull/165) - Fixed consensus variant calling workflow issues (@robert-a-forsyth).
- [#168](https://github.com/IntGenomicsLab/lrsomatic/pull/168) - Fixed default values for `germline_var_keep`, `somatic_var_keep`, `prioritize_caller_germline`, and `prioritize_caller_somatic` parameters to default to `clair` (@robert-a-forsyth).

## v1.0.0 - [28 Nov 2025]

Initial release of IntGenomicsLab/lrsomatic, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- ClairS-TO Module
- cramino module
- modkit pileup module
- mosdepth module
- minimap2/index module
- minimap2/align module
- pigz module
- samtools/cat module
- samtools faidx module
- bam_stats_samtools subworkflow
- mosdepth added to workflow
- add longphase/tag and longphase/phase modules

### `Fixed`

- New channel structure
- No longer possible to have duplicated naming after samtools cat
- restructured minimap2

### `Dependencies`

### `Deprecated`
