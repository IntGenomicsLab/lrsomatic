# IntGenomicsLab/lrsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [2026-04-10]

### `Added`

- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Integrated modkit module for long-read base modification detection and analysis.
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added DeepVariant and DeepSomatic modules for germline and somatic variant calling from long-read sequencing data.
- [#147](https://github.com/IntGenomicsLab/lrsomatic/pull/147) - Implemented whatshap_stats module to generate phase block statistics and phasing quality metrics.
- [#141](https://github.com/IntGenomicsLab/lrsomatic/pull/141) - Added output of phased variants in separate VCF files for improved downstream analysis.
- [#143](https://github.com/IntGenomicsLab/lrsomatic/pull/143) - Added Severus min_support parameter.
- [#145](https://github.com/IntGenomicsLab/lrsomatic/pull/145) - Integrated MultiQC and nanoplot for comprehensive QC reporting with long-read sequencing metrics.
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added GPU support for Clair3, DeepVariant, and fibertools
- [#150](https://github.com/IntGenomicsLab/lrsomatic/pull/150) - Added Claude GitHub Actions workflows for automated code review and PR assistance.

### `Changed`

- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Bulk module versions update. Fixed some issues with Wakhan.
- [#138](https://github.com/IntGenomicsLab/lrsomatic/pull/138) - Perform QC before merging replicates
- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Updated container versions and dependencies for modkit and related tools.
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Refactored variant calling workflow to support both DeepVariant and existing callers with improved configuration handling.
- [#140](https://github.com/IntGenomicsLab/lrsomatic/pull/140) - Improved documentation with additional pipeline usage examples and configuration guidance.


### `Fixed`

- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Resolved Nextflow strict syntax compliance issues for compatibility with latest Nextflow versions.
- [#118](https://github.com/IntGenomicsLab/lrsomatic/pull/118) - Updated nf-core template components to align with latest pipeline standards.
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Corrected bcftools and vcfsplit operations for accurate variant filtering and merging.

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
