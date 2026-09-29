# IntGenomicsLab/lrsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### `Added`

- [#197](https://github.com/IntGenomicsLab/lrsomatic/pull/197) - Added CHM13 support for ClairS-TO's Verdict module, which tags tumour-only calls as germline, somatic or subclonal somatic; its resources were GRCh38-only, so on CHM13 germline variants leaked into `somatic.vcf.gz`. With `--genome CHM13 --skip_ascat` the pipeline builds a CHM13 resource set from the ASCAT files it already downloads and passes it as `--cna_resource_dir`; a prepared directory can be given with `--clairsto_cna_resources` (validated at launch). Without `--skip_ascat` tagging comes from ASCAT's own tables instead (next entry) (@ljwharbers).
- [#197](https://github.com/IntGenomicsLab/lrsomatic/pull/197) - Added `CLAIRSTO_VERDICT_TAG`: when ASCAT is in the run, Verdict's germline tagging is computed from ASCAT's purity, ploidy and segments instead of Verdict's own estimate, so `CLAIRSTO` runs with `--disable_verdict` and ASCAT runs before small variant calling. Output names are unchanged. The tables the tags were computed from are published as `<sample>_Tumor_Purity_Ploidy.txt` and `<sample>_Tumor_CNA.txt`, also on `--skip_ascat` runs (@ljwharbers).
- [#197](https://github.com/IntGenomicsLab/lrsomatic/pull/197) - Added a stub nf-test for `TUMORONLY_SMALLVAR` covering both germline tagging paths (tag `small`) (@ljwharbers).
- [#190](https://github.com/IntGenomicsLab/lrsomatic/pull/190) - Added mutational signature analysis: SigProfilerMatrixGenerator matrices (SBS/DBS/ID) and per-sample COSMIC signature fitting with SigProfilerAssignment on the phased somatic VCF, with CHM13-T2T support via a purpose-built image (`ghcr.io/ljwharbers/sigprofiler`) and CHM13-renormalised COSMIC signatures; new `--sigprofiler_*` / `--download_sigprofiler_genome` / `--skip_signatures` parameters (@ljwharbers).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - Added `modkit_args` (default `--cpg --modified-bases 5mC`) to control the arguments passed to `modkit pileup`, and `modkit_phased` to run the pileup with `--phased` on the Longphase-haplotagged BAMs, producing `_hp1`, `_hp2` and `_combined` bedMethyl files per sample (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added `LRSOMATICREPORT` as the final pipeline step: a self-contained per-sample HTML report covering small variants, structural variants, copy number and QC. Skip it with `--skip_report`; choose the gene panel selected on load with `--report_gene_panel` (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added the [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) tool to the pipeline, first vendored as source at `assets/lrsomatic_report`, now shipped as a container (see the entry below) (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added a `solution_dirs` output to the WAKHAN module so its per-solution copy-number plots can be staged downstream (@ljwharbers).
- [#193](https://github.com/IntGenomicsLab/lrsomatic/pull/193) - Added VEP plugins: AlphaMissense, SIFT/PolyPhen, ClinVar and REVEL on GRCh38, and AlphaMissense plus SIFT/PolyPhen on CHM13 via protein-space lookup. Enabled by default with `--genome GRCh38` or `--genome CHM13` (a first GRCh38 run downloads around 1.4 GB); any `--vep_*` path overrides its default and `--skip_vep_plugins` turns the set off. CADD and EVE are opt-in behind `--vep_cadd_snv`/`--vep_cadd_indel` and `--vep_eve` because of their size (81 GB and 9.6 GB); prepared REVEL and EVE files are published to `<outdir>/vep_plugins/` for reuse. Two lab-hosted AlphaMissense files are CC BY 4.0 with attribution in `CITATIONS.md` (@AmberVerhasselt).
- [#189](https://github.com/IntGenomicsLab/lrsomatic/pull/189) - Added SAVANA structural variant and copy-number calling, running alongside Severus/ASCAT (@yannvrb).

### `Changed`

- [#201](https://github.com/IntGenomicsLab/lrsomatic/pull/201) - `CLAIRSTO` and `CLAIRSTO_VERDICT_TAG` now pull the fork image from Docker Hub: `oras://docker.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13-c0687e8-flat` under Singularity/Apptainer and `docker.io/ljwharbers/clairs-to:0.5.1-verdict-chm13-c0687e8-flat` otherwise, instead of `ghcr.io/ljwharbers/clairs-to`. The `-cpu` SIF on ghcr failed with `PROTOCOL_ERROR` on slow links: ghcr redirects every blob download to an Azure URL that expires at the next 5-minute mark and resets a stream still open then, and Apptainer resumes neither an `oras://` nor a `docker://` download. Docker Hub's download URLs are valid for 50 minutes and only checked when the request starts. `-flat` is the same software copied into an empty image in a few layers (3.3 GB instead of 7 GB); the software and its outputs are unchanged. `docs/usage.md` describes `pullTimeout`, Docker Hub's anonymous pull limit, pre-pulling, and how to recover the remaining `oras://ghcr.io` SIFs resumably (@ljwharbers).
- [#199](https://github.com/IntGenomicsLab/lrsomatic/pull/199) - `CLAIRSTO` and `CLAIRSTO_VERDICT_TAG` now run the `-cpu` rebuild of the fork image (`0.5.1-verdict-chm13-c0687e8-cpu`), which swaps PyTorch's CUDA build for the CPU build of the same version. The software is otherwise unchanged, but the Apptainer SIF drops from 6.53 GB to 3.46 GB. The old image could not be pulled on a normal VSC link: Apptainer fetches an `oras://` SIF as a single unresumable stream, and the signed blob URL ghcr redirects to expires on a 15-minute wall-clock boundary, so 6.53 GB needed 7.3 MB/s sustained and was otherwise cut mid-transfer with `PROTOCOL_ERROR` (@ljwharbers).
- [#197](https://github.com/IntGenomicsLab/lrsomatic/pull/197) - `CLAIRSTO` now runs `ghcr.io/ljwharbers/clairs-to:0.5.1-verdict-chm13-c0687e8` (ClairS-TO 0.5.1) instead of `docker.io/hkubal/clairs-to:v0.4.2`: a fork that lets Verdict read its CNA resources from `--cna_resource_dir`, fixes four places where Verdict's Python port of ASCAT departed from R, and disables Verdict with a warning when its resources cannot be read. **GRCh38 results move as well as CHM13 ones.** Revert to the upstream image once HKU-BAL/ClairS-TO carries these changes. The module also selects the SIF under `-profile apptainer` and sets explicit output prefixes (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - `LRSOMATICREPORT` now runs [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) v1.6.0 from `ghcr.io/ljwharbers/lrsomatic-report:1.6.0` (`oras://ghcr.io/ljwharbers/lrsomatic-report-sif:1.6.0` under Singularity/Apptainer) instead of the source vendored at `assets/lrsomatic_report`, which is deleted; updating the tool is now a tag bump. The image is `linux/amd64` only and **Conda is not supported for this step**: the module errors under `-profile conda`/`mamba`, `conda` is dropped from the nf-test CI matrix, and `--skip_report` gives a conda run (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - Builtin gene panels moved from the tool into this repository at `assets/gene_lists/` and reach the report through `--gene-lists-dir`, so adding a panel needs no report-tool release; `--report_gene_panel` validation reads that directory (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - Removed `--report_src`; render with an unreleased tool by pinning a container tag instead (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - The report's small-variant table now shows the VEP plugin annotations from #193 (AlphaMissense, ClinVar, CADD, REVEL, EVE, and SIFT/PolyPhen on CHM13) as separate class and score columns, present only when the annotated VCF declared the field, with an "Annotation sources" footnote. The v1.4.0 to v1.6.0 releases also fix the SV table's `panel_hit` column, add the optional `applies_to` panel column (next entry) and make the footer read the tool version at render time (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - Rebuilt the builtin `lymphoid` panel from two curated NHL lists (234 rows, up from 72) and scoped its genes per table with the new optional `applies_to` column (`snv`, `sv` or blank for both). Four upstream rows were corrected (`EWSR1`, `KDM6B` and `SIK3` hg38 coordinates; `PRKBC` to `PRKCB`; `RCK` to `DDX6`; duplicate `FAM46C` dropped); see `assets/gene_lists/README.md`. `sarcoma` is untouched (@ljwharbers).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - `MODKIT_PILEUP` now runs `modkit pileup` with `--cpg --modified-bases 5mC` by default; the previous unfiltered pileup produced 30-45 GB per sample, so reruns report only 5mC at CpG sites. Set `--modkit_args=''` (the `=` form; `--modkit_args ''` reaches the pipeline as `true`) or an empty `modkit_args` entry in a params file to get the previous output (@ljwharbers).
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

- [#206](https://github.com/IntGenomicsLab/lrsomatic/pull/206) - A remote (http, https or ftp) ClinVar is now downloaded once per run by the new `VEPPLUGIN_CLINVAR` step in `PREPARE_VEP_PLUGINS`, instead of being staged by Nextflow as a foreign file; local and cloud-storage paths are staged as before. `GERMLINE_VEP` and `SOMATIC_VEP` re-checked the foreign file on its host for every sample, and NCBI answered the burst from a multi-sample GRCh38 run with HTTP 503, so `SOMATIC_VEP` failed with `Can't stage file ...clinvar_20260829.vcf.gz`; `-resume` could not recover, since the failed check changed the staging cache key. The download is checked against the new `--vep_clinvar_md5` and `--vep_clinvar_tbi_md5`, set by default to the checksums NCBI (GRCh38, VCF only) and Ensembl (CHM13, VCF and index) publish, so the pinned release cannot change silently. Resuming a run that already finished re-runs `GERMLINE_VEP` and `SOMATIC_VEP` once, since ClinVar now comes from a task rather than the stage cache. The ClinVar sizes in `docs/usage.md` are also corrected, and `docs/output.md` now documents `vep_plugins/` (@AmberVerhasselt).
- [#203](https://github.com/IntGenomicsLab/lrsomatic/pull/203) - `CLAIRS` no longer runs with `--haplotagged_tumor_bam_provided_so_skip_intermediate_phasing_and_haplotagging`. Since somatic calling was moved ahead of `PHASING_HAPLOTYPING` (v1.1.0), ClairS has received the untagged minimap2 BAM, so the flag made it skip its own phasing and haplotagging and call every paired sample without haplotype information: the full-alignment model saw no `HP` tags and the haplotype filtering step had nothing to filter on, the same as `--disable_phasing`. ClairS now runs Clair3 on the normal and tumour BAMs and phases the tumour itself. **Paired somatic calls change** (fewer false positives expected), and `CLAIRS` takes longer and uses more work-directory space (@ljwharbers).
- [#203](https://github.com/IntGenomicsLab/lrsomatic/pull/203) - `docs/output.md` now lists the ClairS SNV output as `snvs.vcf.gz`, the name the pipeline publishes, instead of `snv.vcf.gz` (@ljwharbers).
- [#196](https://github.com/IntGenomicsLab/lrsomatic/pull/196) - `LRSOMATICREPORT` now points `XDG_CACHE_HOME` at the task directory alongside `HOME` and `TMPDIR`. Singularity/Apptainer inherit the host environment, so on sites that set it outside the bind-mounted work tree the render died with `Read-only file system (os error 30): mkdir '<...>/.cache/quarto'` (@AmberVerhasselt, @ljwharbers).
- [#193](https://github.com/IntGenomicsLab/lrsomatic/pull/193) - `--vep_eve https://evemodel.org/api/proteins/bulk/download/` was rejected at launch because the "needs preparing" check keyed on a `.zip` suffix; it now checks whether the value is already a prepared bgzipped file (@AmberVerhasselt).
- [#188](https://github.com/IntGenomicsLab/lrsomatic/pull/188) - `MODKIT_PILEUP` now runs a patched modkit 0.6.4 ([ljwharbers/modkit@pacbio-conflict-fix](https://github.com/ljwharbers/modkit/tree/pacbio-conflict-fix)): `ghcr.io/ljwharbers/modkit:0.6.4-pacbiofix-6e0afa2` under Docker and `oras://ghcr.io/ljwharbers/modkit-sif:0.6.4-pacbiofix-6e0afa2` under Singularity/Apptainer. Stock modkit 0.4.3-0.6.4 dropped 32-65 % of reads from recent PacBio HiFi BAMs and returned empty `--cpg` pileups ([nanoporetech/modkit#612](https://github.com/nanoporetech/modkit/issues/612), fix proposed in [nanoporetech/modkit#720](https://github.com/nanoporetech/modkit/pull/720)), and ignored `--phased`/`--modified-bases` for PacBio BAMs with 6mA calls. The image is `linux/amd64` only and Conda is not supported (use `--skip_modkit` there); return to the biocontainer once a release includes the fix (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - NanoPlot pre- and post-alignment statistics now reach MultiQC; their channels were declared but never assigned, so the NanoStat section was silently missing from every report (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - `NANOPLOT_PRE` now includes the replicate in its output prefix, so samples with more than one replicate no longer collapse into a single MultiQC sample (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - samtools stats/flagstat/idxstats and mosdepth outputs are now prefixed `<sample>_<type>` instead of `<sample>`; the tumor and normal files of a matched pair shared a basename, so MultiQC logged `Duplicate sample name found! Overwriting` and reported only one of them (@ljwharbers).
- [#194](https://github.com/IntGenomicsLab/lrsomatic/pull/194) - A sample whose `clair3_model` differs from the model implied by its BAM header no longer makes Clair3 run twice on every normal BAM sharing that header model, with a random copy winning the downstream join. The same bug made `UNTAR` fail with `mkdir: missing operand` when the header's basecall model was absent from the model map ([#191](https://github.com/IntGenomicsLab/lrsomatic/pull/191)). Model entries are now keyed on the downloaded model (@YannVRB, @ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - Raised MultiQC's `log_filesize_limit` to 500 MB, because long-read `samtools stats` output exceeds the 50 MB default and was silently skipped (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - The `clair_only` nf-test now asserts sample4's merged tumour BAM with `bam().getReadsMD5()` instead of file md5s, which `samtools merge` makes irreproducible by giving colliding `@PG` IDs a random suffix (@ljwharbers).
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
