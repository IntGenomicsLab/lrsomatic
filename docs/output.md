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
├── pipeline_info
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

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample.before_correction.sample.tumour.germline.png`  | LogR and BAF plots from the normal sample |
| `sample.before_correction.sample.tumour.tumour.png` | LogR and BAF plots from the tumor sample
| `sample.cnvs.txt` | a tsv file describing each chromosome segment with a copy number alteration and it's major and minor copy number
| `sample.metrics.txt` | a tsv file describing summary statistics for the sample 
| `sample.normal_alleleFrequencies_chr(1-22,X).txt` | a tsv file describing the snp counts for the normal sample at each position and their respective depths
| `sample.purityploidy.txt` | a tsv file describing the purity and ploidy values of the sample
| `sample.segments.txt` | a tsv file describing each chromosome segment and it's major and minor copy number
| `sample.tumour_alleleFrequencies_chr(1-22,X).txt` | a tsv file describing the snp counts for the tumor sample at each position and their respective depths
| `sample.tumour_normalBAF_rawBAF.txt` | a tsv file with the raw BAF values in the normal sample
| `sample.tumour_normalBAF.txt` | a tsv file with the BAF values in the normal sample
| `sample.tumour_normalLogR.txt` | a tsv file with the LogR values in the normal sample
| `sample.tumour_tumourBAF_rawBAF.txt` | a tsv file with the raw BAF values in the tumor sample
| `sample.tumour_tumourBAF.txt` | a tsv file with the BAF values in the tumor sample
| `sample.tumour_tumourLogR.txt` | a tsv file with the LogR values in the tumor sample
| `sample.tumour.ASCATprofile.png` | a png file with the corrected overall copy number profile with ploidy, purity, and goodness of fit metrics
| `sample.tumour.ASPCF.png` | a png file with the correct LogR and BAF plots of the tumor sample
| `sample.tumour.rawprofile.png` | a png file with the raw overall copy number profile with ploidy, purity, and goodness of fit metrics
| `sample.tumour.sunrise.png` | a png file with a purity and ploidy fit

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

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`sample_normal.bam` | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the normal sample
| `sample_normal.bam.bai` | index file for the normal bam file
| ` sample_tumor.bam` | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the tumor sample
| `sample_tumor.bam.bai` | index file for the tumor bam file
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

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`cramino_aln/sample.cramino.txt` | cramino QC summary statistics for the aligned bam file
| `cramino_ubam/sample.cramino.txt` | cramino QC summary statistics for the unaligned bam files
| `fibertoolsrs/sample_qc.txt` | fibertools QC summary for the bam file
| `mosdepth/sample.mosdepth.global.dist.txt` | a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value
| `mosdepth/sample.mosdepth.summary.txt` | overall summary file from mosdepth tool
|`samtools/sample.flagstat` |
|`samtools/sample.idxstats` |
|`samtools/sample.stats` |
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
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`merge_output.vcf.gz` | Merged germline indel and snv calls in vcf format
|`merge_output.vcf.gz` | index for germline small variant calls

#### `clairS`
```
├── clairS
│   ├── indel.vcf.gz
│   ├── indel.vcf.gz.tbi
│   ├── snv.vcf.gz
│   ├── snv.vcf.gz.tbi
```
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`indel.vcf.gz` | Somatic indel calls in vcf format
|`indel.vcf.gz.tbi` | Index for somatic indel calls
|`snv.vcf.gz` | Somatic SNV calls in vcf format
|`snv.vcf.gz.tbi`| Index for somatic SNV calls

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

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`germline.vcf.gz` | SNV and indel calls marked as germline (will not include variants QC)
|`germline.vcf.gz.tbi` | Index file for germline small variant calls
|`indel.vcf.gz` | Raw indel calls in vcf format
|`indel.vcf.gz.tbi` | Index for somatic indel calls
|`snv.vcf.gz` | Raw SNV calls in vcf format
|`snv.vcf.gz.tbi`| Index for  SNV calls
|`somatic.vcf.gz`| SNV and indel calls marked as PASS and without a germline tag
|`somatic.vcf.gz`| Index for osmatic small variatn calls

#### `severus`
```
├── severus
│   ├── all_SVs
│   │   ├── plots
│   │   │   ├──
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_all.vcf.gz
│   │   ├── severus_all.vcf.gz.tbi
│   ├── somatic_SVs
│   │   ├── plots
│   │   │   ├──
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_somatic.vcf.gz
│   │   ├── severus_somatic.vcf.gz.tbi
│   ├── breakpoints_double.csv
│   ├── read_ids.csv
│   ├── read_qual.txt
│   ├── severus.log
```
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`all_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in all clustered events
|`all_SVs/breakpoint_cluster.tsv` | a tsv containing all clustered events
|`all_SVs/severus_all.vcf.gz` | A vcf file containing all identified structural variants
|`all_SVs/severus_all.vcf.gz.tbi` | Index for all identified structural variants
|`somatic_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in somatic clustered events
|`somatic_SVs/breakpoint_cluster.tsv` | a tsv containing somatic clustered events
|`somatic_SVs/severus_somatic.vcf.gz` | A vcf file containing identified somatic structural variants
|`somatic_SVs/severus_somatic.vcf.gz.tbi` | Index for identified somatic structural variants



</details>

### `vep`

### `pipeline_info`
## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [samtools/cat](#samtools_cat) - Concatenate replicates from the same sample
- [Cramino](#cramino) - Quality check the aligned or unaligned lomg read bam files
- [fibertools/predictm6a](#fibertools) -
- [fibertools/nucleosome](#fibertools)-
- [fibertools/fire](#fibertools) - 
- [fibertools/qc](#fibertools) -
- [minimap2/index](#minimap2) - Creating minimap2 index
- [minimap2/align](#minimap2) - Aligning long-read samples to the reference genome
- [clair3](#clair3) - 
- [clairsto](#clairsto)-
- [clairs](#clairs) - 
- [longphase/phase](#longphase) -
- [longphase/haplotag](#longphase) -
- [severus](#severus) - 
- [wakhan](#wakhan) - 
- [Modkit](#modkit) - Tool to process methylation data
- [Mosdepth](#mosdepth) - Tool to check the sequencing depth
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### samtools

### Minimap2

[ClairS-TO](https://github.com/HKU-BAL/ClairS-TO) ool in the Clair series to support long-read somatic small variant calling with only tumor samples available. For more information, see <https://github.com/HKU-BAL/ClairS-TO>.

### Cramino

<details markdown="1">
<summary>Output files</summary>

- `cramino/`
  - `*_cramino.txt`: text file containing the quality check results from Cramino

</details>

[Cramino](https://github.com/wdecoster/cramino) is a tool for fast quality check of aligned or unaligned long read sequencing bam files. For more information, check <https://github.com/wdecoster/cramino>.

### Modkit

<details markdown="1">
<summary>Output files</summary>

- `modkit/`
  - `*.bed`: bed file with genomic positions with modified bases
  - `*.bedgraph`: bedgraph file with genomic positions with modified bases
  - `*log`: log file

</details>

[Modkit](https://github.com/nanoporetech/modkit) is a tool to work with methylated bases from bam files. It allows you to update the bam file modification information, filter the information, or to extract it to other file formarts such as bed file. For more information, see <https://github.com/nanoporetech/modkit>.

### Mosdepth

<details markdown="1">
<summary>Output files</summary>

- `mosdepth`
  - `*.summary.txt'`: summary text file
  - `*.global.dist.txt`: global information on sequencing depth file
  - optionally outputs other files by regions.

</details>

[Mosdepth](https://github.com/brentp/mosdepth) is a tool to work with WGS/exome/targeted sequencing data to obtain sequecing depth information, For more information, see <https://github.com/brentp/mosdepth>.

### ClairS-TO

<details markdown="1">
<summary>Output files</summary>

- `clairsto/`
  - `germline.vcf.gz`: a VCF file containing germline SNP and indel calls
  - `somatic.vcf.gz`: a VCF file containing somatic SNP and indel calls
  - `logs`: directory containing log files for the ClairS-TO run
  - `run_clairs_to.log`: a text file containing a summary of the ClairS-TO run
  - `run_clairs_to.log.bak`: a text file containing the run command

</details>

### Severus

<details markdown="1">
<summary>Output files</summary>

- `severus`
- `severus.log`: severus log file
- `read_qual.txt`: the read quality
- `breakpoints_double.csv`: the breakpoint file

</details>

[Severus](https://github.com/KolmogorovLab/Severus) is a tool to call not only somatic, but also germline structural variant calls. For mor information on the tool and its usage, check out <https://github.com/KolmogorovLab/Severus>.

### Longphase

<details markdown="1">
<summary>Output files</summary>

TODO - add description ot output file

- `longphase/`

</details>

[Longphase](https://github.com/twolinin/longphase) isa tool to phase your input variant calls and haplotag your bam. To see more, see <https://github.com/twolinin/longphase>.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
