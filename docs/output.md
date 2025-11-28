# IntGenomicsLab/lrsomatic: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

### Output Example

```
├── Sample 1
│    ├── ascat
│    ├── bamfiles
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    ├── variants
│    │   ├──clairS-TO
│    │   ├──severus
│    ├── vep
│    │   ├── germline
│    │   ├── somatic
│    │   ├── SVs
│
├── Sample 2
│    ├── ascat
│    ├── bamfiles
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    │    ├── normal
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    ├── variants
│    │   ├── clair3
│    │   ├── clairS
│    │   ├── severus
│    ├── vep
│    │   ├── germline
│    │   ├── somatic
│    │   ├── SVs
│    ├── wakhan
├── pipeline_info
├── multiqc

```

### `ascat`

<details markdown="1">
<summary>Output files</summary>

```
├── ascat
│   ├── sample.before_correction.sample.tumour.germline.png
│   ├── sample.before_correction.sample.tumour.tumour.png
│   ├── sample.cnvs.txt
│   ├── sample.metrics.txt
│   ├── sample.normal_alleleFrequencies_chr(1-22,X).txt
│   ├── sample.purityploidy.txt
│   ├── sample.segments.txt
│   ├── sample.tumour_alleleFrequencies_chr(1-22,X).txt
│   ├── sample.tumour_normalBAF_rawBAF.txt
│   ├── sample.tumour_normalBAF.txt
│   ├── sample.tumor_tumourLogR.txt
│   ├── sample.tumour.ASCATprofile.png
│   ├── sample.tumour.ASPCF.png
│   ├── sample.tumour.rawprofile.png
│   ├── sample.tumour.sunrise.png
```

| File                                                  | Description                                                                                                      |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `sample.before_correction.sample.tumour.germline.png` | LogR and BAF plots from the normal sample before correction                                                      |
| `sample.before_correction.sample.tumour.tumour.png`   | LogR and BAF plots from the tumor sample before correction                                                       |
| `sample.cnvs.txt`                                     | a tsv file describing each chromosome segment with a copy number alteration and it's major and minor copy number |
| `sample.metrics.txt`                                  | a tsv file describing summary statistics for the sample                                                          |
| `sample.normal_alleleFrequencies_chr(1-22,X).txt`     | a tsv file describing the snp counts for the normal sample at each position and their respective depths          |
| `sample.purityploidy.txt`                             | a tsv file describing the purity and ploidy values of the sample                                                 |
| `sample.segments.txt`                                 | a tsv file describing each chromosome segment and it's major and minor copy number                               |
| `sample.tumour_alleleFrequencies_chr(1-22,X).txt`     | a tsv file describing the snp counts for the tumor sample at each position and their respective depths           |
| `sample.tumour_normalBAF_rawBAF.txt`                  | a tsv file with the raw BAF values in the normal sample                                                          |
| `sample.tumour_normalBAF.txt`                         | a tsv file with the BAF values in the normal sample                                                              |
| `sample.tumour_normalLogR.txt`                        | a tsv file with the LogR values in the normal sample                                                             |
| `sample.tumour_tumourBAF_rawBAF.txt`                  | a tsv file with the raw BAF values in the tumor sample                                                           |
| `sample.tumour_tumourBAF.txt`                         | a tsv file with the corrected BAF values in the tumor sample                                                     |
| `sample.tumour_tumourLogR.txt`                        | a tsv file with the corrected LogR values in the tumor sample                                                    |
| `sample.tumour.ASCATprofile.png`                      | a png file with the corrected overall copy number profile with ploidy, purity, and goodness of fit metrics       |
| `sample.tumour.ASPCF.png`                             | a png file with the corrected LogR and BAF plots of the tumor sample                                             |
| `sample.tumour.rawprofile.png`                        | a png file with the raw overall copy number profile with ploidy, purity, and goodness of fit metrics             |
| `sample.tumour.sunrise.png`                           | a png file with a purity and ploidy fit                                                                          |

</details>

### `bamfiles`

<details markdown="1">
<summary>Output files</summary>

```
├── bamfiles
│   ├── sample_normal.bam
│   ├── sample_normal.bam.bai
│   ├── sample_tumor.bam
│   ├── sample_tumor.bam.bai
```

| File                    | Description                                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------------------- |
| `sample_normal.bam`     | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the normal sample |
| `sample_normal.bam.bai` | index file for the normal bam file                                                                   |
| ` sample_tumor.bam`     | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the tumor sample  |
| `sample_tumor.bam.bai`  | index file for the tumor bam file                                                                    |

</details>

### `qc`

<details markdown="1">
<summary>Output files</summary>

```
├── qc
│   ├── tumor
│   │   ├── cramino_aln
│   │   │   ├── sample.cramino.txt
│   │   ├── cramino_ubam
│   │   │   ├── sample.cramino.txt
│   │   ├── fibertoolsrs
│   │   │   ├── sample_qc.txt
│   │   ├── mosdepth
│   │   │   ├── sample.mosdepth.global.dist.txt
│   │   │   ├── sample.mosdepth.summary.txt
│   │   ├── samtools
│   │   │   ├── sample.flagstat
│   │   │   ├── sample.idxstats
│   │   │   ├── sample.stats
```

| File                                       | Description                                                                                                              |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `cramino_aln/sample.cramino.txt`           | cramino QC summary statistics for the aligned bam file                                                                   |
| `cramino_ubam/sample.cramino.txt`          | cramino QC summary statistics for the unaligned bam files                                                                |
| `fibertoolsrs/sample_qc.txt`               | fibertools QC summary for the bam file                                                                                   |
| `mosdepth/sample.mosdepth.global.dist.txt` | a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value |
| `mosdepth/sample.mosdepth.summary.txt`     | overall summary file from mosdepth tool                                                                                  |
| `samtools/sample.flagstat`                 | a summary of the counts of different samtools flags                                                                      |
| `samtools/sample.idxstats`                 | a summary of the number of mapped and unmapped reads                                                                     |
| `samtools/sample.stats`                    | summary statistics from the bamfile                                                                                      |

</details>

### `variants`

<details markdown="1">
<summary>Output files</summary>

#### `clair3`

```
├── clair3
│   ├── merge_output.vcf.gz
│   ├── merge_output.vcf.gz.tbi
```

| File                  | Description                                       |
| --------------------- | ------------------------------------------------- |
| `merge_output.vcf.gz` | Merged germline indel and snv calls in vcf format |
| `merge_output.vcf.gz` | index for germline small variant calls            |

#### `clairS`

```
├── clairS
│   ├── indel.vcf.gz
│   ├── indel.vcf.gz.tbi
│   ├── snv.vcf.gz
│   ├── snv.vcf.gz.tbi
```

| File               | Description                       |
| ------------------ | --------------------------------- |
| `indel.vcf.gz`     | Somatic indel calls in vcf format |
| `indel.vcf.gz.tbi` | Index for somatic indel calls     |
| `snv.vcf.gz`       | Somatic SNV calls in vcf format   |
| `snv.vcf.gz.tbi`   | Index for somatic SNV calls       |

#### `clairS-TO`

```
├── clairS
│   ├── germline.vcf.gz
│   ├── germline.vcf.gz.tbi
│   ├── indel.vcf.gz
│   ├── indel.vcf.gz.tbi
│   ├── snv.vcf.gz
│   ├── snv.vcf.gz.tbi
│   ├── somatic.vcf.gz
│   ├── somatic.vcf.gz.tbi
```

| File                  | Description                                                           |
| --------------------- | --------------------------------------------------------------------- |
| `germline.vcf.gz`     | SNV and indel calls marked as germline (will not include variants QC) |
| `germline.vcf.gz.tbi` | Index file for germline small variant calls                           |
| `indel.vcf.gz`        | Raw indel calls in vcf format                                         |
| `indel.vcf.gz.tbi`    | Index for somatic indel calls                                         |
| `snv.vcf.gz`          | Raw SNV calls in vcf format                                           |
| `snv.vcf.gz.tbi`      | Index for SNV calls                                                   |
| `somatic.vcf.gz`      | SNV and indel calls marked as PASS and without a germline tag         |
| `somatic.vcf.gz`      | Index for osmatic small variatn calls                                 |

#### `severus`

```
├── severus
│   ├── all_SVs
│   │   ├── plots
│   │   │   ├── severus_{*}.html
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_all.vcf.gz
│   │   ├── severus_all.vcf.gz.tbi
│   ├── somatic_SVs
│   │   ├── plots
│   │   │   ├── severus_{*}.html
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_somatic.vcf.gz
│   │   ├── severus_somatic.vcf.gz.tbi
│   ├── breakpoints_double.csv
│   ├── read_ids.csv
│   ├── read_qual.txt
│   ├── severus.log
```

| File                                      | Description                                                                       |
| ----------------------------------------- | --------------------------------------------------------------------------------- |
| `all_SVs/plots/severus_{*}.html`          | html file containing a plot of connected breakpoints in a cluster                 |
| `all_SVs/breakpoint_cluster_list.tsv`     | tsv containing the breakpoints in all clustered events                            |
| `all_SVs/breakpoint_cluster.tsv`          | a tsv containing all clustered events                                             |
| `all_SVs/severus_all.vcf.gz`              | A vcf file containing all identified structural variants                          |
| `somatic_SVs/plots/severus_{*}.html`      | html file containing a plot of connected breakpoints in a cluster                 |
| `somatic_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in somatic clustered events                        |
| `somatic_SVs/breakpoint_cluster.tsv`      | a tsv containing somatic clustered events                                         |
| `somatic_SVs/severus_somatic.vcf.gz`      | A vcf file containing identified somatic structural variants                      |
| `somatic_SVs/severus_somatic.vcf.gz.tbi`  | Index for identified somatic structural variants                                  |
| `breakpoints_double.csv`                  | csv file containing detailed information about identified breakpoints in bam file |
| `read_ids.csv`                            | a csv file containing read ids associated with each identified SV                 |
| `read_qual.txt`                           | file containing quality statistics about identified segements                     |
| `severus.log`                             | log file                                                                          |

</details>

### `vep`

<details markdown="1">
<summary>Output files</summary>

```
├── vep
│   ├── germline
│   │   ├── sample_GERMLINE_VEP.vcf.gz
│   │   ├── sample_GERMLINE_VEP_summary.html
│   │   ├── sample_GERMLINE_VEP.vcf.gz.tbi
│   ├── somatic
│   │   ├── sample_SOMATIC_VEP.vcf.gz
│   │   ├── sample_SOMATIC_VEP_summary.html
│   │   ├── sample_SOMATIC_VEP.vcf.gz.tbi
│   ├── SVs
│   │   ├── sample_SV_VEP.vcf.gz
│   │   ├── sample_SV_VEP_summary.html
│   │   ├── sample_SV_VEP.vcf.gz.tbi
```

| File                                        | Description                                                             |
| ------------------------------------------- | ----------------------------------------------------------------------- |
| `germline/sample_GERMLINE_VEP.vcf.gz`       | Annotated germline indel and SNV vcf file                               |
| `germline/sample_GERMLINE_VEP_summary.html` | Visual summary of germline indel and SNV annotations in html format     |
| `germline/sample_GERMLINE_VEP.vcf.gz.tbi`   | Annotated germline indel and SNV vcf index file                         |
| `somatic/sample_SOMATIC_VEP.vcf.gz`         | Annotated somatic indel and SNV vcf file                                |
| `somatic/sample_SOMATIC_VEP_summary.html`   | Visual summary of somatic indel and SNV annotations in html format      |
| `somatic/sample_SOMATIC_VEP.vcf.gz.tbi`     | Annotated somatic indel and SNV vcf index file                          |
| `SVs/sample_SV_VEP.vcf.gz`                  | Annotated somatic structural variant vcf file                           |
| `SVs/sample_SV_VEP_summary.html`            | Visual summary of somatic structural variant annotations in html format |
| `SVs/sample_SV_VEP.vcf.gz.tbi`              | Annotated somatic structural variant vcf index file                     |

</details>

### `wakhan`

<details markdown="1">
<summary>Output files</summary>

```
├── wakhan
│   ├── {ploidy}_{purity}_{confidence}
│   │   ├── bed_output
│   │   │   ├── genes_copynumber_states.bed
│   │   │   ├── loh_regions.bed
│   │   │   ├── sample_{ploidy}_{purity}_{confidence}_HP_1.bed
│   │   │   ├── sample_{ploidy}_{purity}_{confidence}_HP_2.bed
│   │   ├── variation_plots
│   │   │   ├── chr{1-22,X,Y}_cn.html
│   │   │   ├── chr{1-22,X,Y}_cn.pdf
│   │   │   ├── CN_VARIATION_INDEX.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genes_genome.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genes_genome.pdf
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_breakpoints.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_breakpoints.pdf
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.pdf
│   ├── coverage_data
│   │   ├── {0-23}_SNPS.csv
│   │   ├── coverage_ps.csv
│   │   ├── phase_corrected_coverage.csv
│   │   ├── pileup_SNPs.csv
│   ├── coverage_plots
│   │   ├── chr{1-22,X,Y}_cov.html
│   │   ├── chr{1-22,X,Y}.pdf
│   │   ├── COVERAGE_INDEX.html
│   ├── phasing output
│   │   ├── chr{1-22,X,Y}_phase_correction_0.html
│   │   ├── chr{1-22,X,Y}_phase_correction_1.html
│   │   ├── chr{1-22,X,Y}_without_phase_correction.html
│   │   ├── chr{1-22,X,Y}.pdf
│   │   ├── sample.rephased.vcf.gz
│   │   ├── sample.rephased.vcf.gz.tbi
│   ├── sample_heatmap_ploidy_purity.html
│   ├── sample_heatmap_ploidy_purity.html.pdf
│   ├── sample_optimized_peak.html
│   ├── solutions_ranks.tsv

```

| File                                                                                                   | Description                                                                          |
| ------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ |
| `{ploidy}_{purity}_{confidence}/bed_output/genes_copynumber_states.bed`                                | bed file containing allele specific copy number values with coverage information     |
| `{ploidy}_{purity}_{confidence}/bed_output/loh_regions.bed`                                            | bed file containing positions of loss of heterozygosity regions                      |
| `{ploidy}_{purity}_{confidence}/bed_output/sample_{ploidy}_{purity}_{confidence}_HP_1.bed`             | bed file containing copy number states, coverage, and SV breakpoints for haplotype 1 |
| `{ploidy}_{purity}_{confidence}/bed_output/sample_{ploidy}_{purity}_{confidence}_HP_2.bed`             | bed file containing copy number states, coverage, and SV breakpoints for haplotype 2 |
| `{ploidy}_{purity}_{confidence}/variation_plots/chr{1-22,X,Y}_cn.html`                                 | html based plotly plot of copy number and coverage for individual chromosomes        |
| `{ploidy}_{purity}_{confidence}/variation_plots/chr{1-22,X,Y}_cn.pdf`                                  | pdf based plotly plot of copy number and coverage for individual chromosomes         |
| `{ploidy}_{purity}_{confidence}/variation_plots/CN_VARIATION_INDEX.html`                               | unclear html plot                                                                    |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genes_genome.html`               | html plots of copy number variations in highlighted genes                            |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genes_genome.pdf`                | pdf plots of copy number variations in highlighted genes                             |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.html` | genome-wide html copy number plots with coverage information on same axis            |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.pdf`  | genome-wide pdf copy number plots with coverage information on same axis             |
| `coverage_data/{0-23}_SNP.csv`                                                                         | CSV of coverage data per chromosome                                                  |
| `coverage_data/coverage_ps.csv`                                                                        | CSV of overall haplotype specific coverage data                                      |
| `coverage_data/coverage.csv`                                                                           | CSV of overall coverage data                                                         |
| `coverage_data/phase_corrected_coverage.csv`                                                           | CSV of overall phase-corrected coverage data                                         |
| `coverage_data/pileup_SNPs.csv`                                                                        | CSV of SNP pileup data                                                               |
| `coverage_plots/chr{1-22,X,Y}_cov.html`                                                                | chromosome specific html coverage plots                                              |
| `coverage_plots/chr{1-22,X,Y}_cov.pdf`                                                                 | chromosome specific pdf coverage plots                                               |
| `coverage_plots/COVERAGE_INDEX.html`                                                                   | unclear html plot                                                                    |
| `phasing_output/chr{1-23,X,Y}_phase_correction_0.html`                                                 | Phase-switch error correction plot per chromosome                                    |
| `phasing_output/chr{1-23,X,Y}_phase_correction_1.html`                                                 | Phase-switch error correction plot per chromosome                                    |
| `phasing_output/chr{1-22,X,Y}_without_phase_correction.html`                                           | Phase-switch error without phase correction plot per chromosome                      |
| `phasing_output/chr{1-22,X,Y}.pdf`                                                                     | Phase-switch error correction plot                                                   |
| `phasing_output/PHASE_CORRECTION_INDEX`                                                                | unclear html plot                                                                    |
| `phasing_output/sample_rephased.vcf.gz`                                                                | phase corrected SNP vcf file                                                         |
| `phasing_output/sample_rephased.vcf.gz.tbi`                                                            | phase corrected SNP vcf index file                                                   |
| `sample_heatmap_ploidy_purity.html`                                                                    | heatmap html plot of purity ploidy fit                                               |
| `sample_heatmap_ploidy_purity.html.pdf`                                                                | heatmap html plot of purity ploidy fit                                               |
| `sample_optimized_peak.html`                                                                           | optimization peak plot                                                               |
| `solutions_ranks.tsv`                                                                                  | rank of potential purity ploidy solutions                                            |

</details>

### `multiqc`

<details markdown="1">
<summary>Output files</summary>

```
├── multiqc
│   ├── multiqc_data
│   │   ├── BETA-multiqc.parquet
│   │   ├── llms-full.txt
│   │   ├── mosdepth_cov_dist.txt
│   │   ├── mosdepth_cumcov_dist.txt
│   │   ├── mosdepth_perchrom.txt
│   │   ├── mosdepth-coverage-per-contig-multi.txt
│   │   ├── mosdepth-cumcoverage-dist-id.txt
│   │   ├── mosdepth-xy-coverage-plot.txt
│   │   ├── multiqc_citations
│   │   ├── multiqc_data.json
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc_software_versions.txt
│   │   ├── multiqc_sources.txt
│   │   ├── multiqc.log
│   ├── multiqc_plots
│   │   ├── pdf
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.pdf
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.pdf
│   │   │   ├── mosdepth-cumcoverage-dist-id.pdf
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.pdf
│   │   │   ├── mosdepth-xy-coverage-plot-pct.pdf
│   │   ├── png
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.png
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.png
│   │   │   ├── mosdepth-cumcoverage-dist-id.png
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.png
│   │   │   ├── mosdepth-xy-coverage-plot-pct.png
│   │   ├── svg
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.svg
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.svg
│   │   │   ├── mosdepth-cumcoverage-dist-id.svg
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.svg
│   │   │   ├── mosdepth-xy-coverage-plot-pct.svg
│   ├── multiqc_report.html

```

| File                                                           | Description                                                  |
| -------------------------------------------------------------- | ------------------------------------------------------------ |
| `multiqc_data/BETA-multiqc.parquet`                            | Multiqc data in Apache Parquet format (BETA)                 |
| `multiqc_data/llms-full.txt`                                   | Prompt for large-language-model summary                      |
| `multiqc_data/mosdepth_cov_dist.txt`                           | text file of coverage distribution                           |
| `multiqc_data/mosdepth_cumcov_dist.txt`                        | text file of cummulative coverage distribution               |
| `multiqc_data/mosdepth_perchrom.txt`                           | text file of coverage per chromosome                         |
| `multiqc_data/mosdepth-coverage-per-contig-multi.txt`          | text file of coverage per contig                             |
| `multiqc_data/mosdepth-cumcoverage-dist-id.txt`                | unclear text file                                            |
| `multiqc_data/mosdepth-xy-coverage-plot.txt`                   | summary of chr X, chr Y coverage                             |
| `multiqc_data/multiqc_citations`                               | citations for multiqc                                        |
| `multiqc_data/multiqc_data.json`                               | json file containing multiqc data output                     |
| `multiqc_data/multiqc_general_stats.txt`                       | summary statistics for the samples                           |
| `multiqc_data/multiqc_software_versions.txt`                   | software versions for tools used in multiqc                  |
| `multiqc_data/multiqc_sources.txt`                             | software sources for tools used in multiqc                   |
| `multiqc_data/multiqc.log`                                     | log file for multiqc                                         |
| `multiqc_plots/pdf/mosdepth-coverage-per-contig-multi-cnt.pdf` | pdf format plot of coverage per contig                       |
| `multiqc_plots/pdf/mosdepth-coverage-per-contig-multi-log.pdf` | pdf format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/pdf/mosdepth-cumcoverage-dist-id.pdf`           | pdf format plot of distribution of cumulative coverage       |
| `multiqc_plots/pdf/mosdepth-xy-coverage-plot-cnt.pdf`          | pdf format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/pdf/mosdepth-xy-coverage-plot-pct.pdf`          | pdf format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/png/mosdepth-coverage-per-contig-multi-cnt.png` | png format plot of coverage per contig                       |
| `multiqc_plots/png/mosdepth-coverage-per-contig-multi-log.png` | png format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/png/mosdepth-cumcoverage-dist-id.png`           | png format plot of distribution of cumulative coverage       |
| `multiqc_plots/png/mosdepth-xy-coverage-plot-cnt.png`          | png format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/png/mosdepth-xy-coverage-plot-pct.png`          | png format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/svg/mosdepth-coverage-per-contig-multi-cnt.svg` | svg format plot of coverage per contig                       |
| `multiqc_plots/svg/mosdepth-coverage-per-contig-multi-log.svg` | svg format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/svg/mosdepth-cumcoverage-dist-id.svg`           | svg format plot of distribution of cumulative coverage       |
| `multiqc_plots/svg/mosdepth-xy-coverage-plot-cnt.svg`          | svg format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/svg/mosdepth-xy-coverage-plot-pct.svg`          | svg format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/multiqc_report.html`                            | svg format plot of chr X and chr Y coverage by percentage    |

### `pipeline_info`

<details markdown="1">
<summary>Output files</summary>

```
├── pipeline_info
│   ├── execution_report_{DATE}.html
│   ├── execution_timeline_{DATE}.html
│   ├── execution_trace_{DATE}.txt
│   ├── lrsomatic_softwar_mqc_versions.yml
│   ├── params_{DATE}.json
│   ├── pipeline_dag_{DATE}/html
```

| File                                 | Description                                                                                 |
| ------------------------------------ | ------------------------------------------------------------------------------------------- |
| `execution_report_{DATE}.hmtl`       | summary of pipeline resource and timing usage in a html report                              |
| `execution_timeline_{DATE}.hmtl`     | a graphical summary of the timing of each module's task over the course of the pipeline run |
| `lrsomatic_softwar_mqc_versions.yml` | summary of the versions of each tool used by the pipeline                                   |
| `params_{DATE}.json`                 | summary of the paramaters used in the pipeline                                              |
| `pipeline_dag_{DATE}.html`           | flow chart summarizing the pipeline run                                                     |

</details>
