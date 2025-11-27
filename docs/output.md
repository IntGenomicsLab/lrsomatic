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

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample.before_correction.sample.tumour.germline.png`  | LogR and BAF plots from the normal sample before correction |
| `sample.before_correction.sample.tumour.tumour.png` | LogR and BAF plots from the tumor sample before correction
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
| `sample.tumour_tumourBAF.txt` | a tsv file with the corrected BAF values in the tumor sample
| `sample.tumour_tumourLogR.txt` | a tsv file with the corrected LogR values in the tumor sample
| `sample.tumour.ASCATprofile.png` | a png file with the corrected overall copy number profile with ploidy, purity, and goodness of fit metrics
| `sample.tumour.ASPCF.png` | a png file with the corrected LogR and BAF plots of the tumor sample
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
|`samtools/sample.flagstat` | a summary of the counts of different samtools flags
|`samtools/sample.idxstats` | a summary of the number of mapped and unmapped reads
|`samtools/sample.stats` | summary statistics from the bamfile
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
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`all_SVs/plots/severus_{*}.html` | html file containing a plot of connected breakpoints in a cluster
|`all_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in all clustered events
|`all_SVs/breakpoint_cluster.tsv` | a tsv containing all clustered events
|`all_SVs/severus_all.vcf.gz` | A vcf file containing all identified structural variants
|`somatic_SVs/plots/severus_{*}.html` | html file containing a plot of connected breakpoints in a cluster
|`somatic_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in somatic clustered events
|`somatic_SVs/breakpoint_cluster.tsv` | a tsv containing somatic clustered events
|`somatic_SVs/severus_somatic.vcf.gz` | A vcf file containing identified somatic structural variants
|`somatic_SVs/severus_somatic.vcf.gz.tbi` | Index for identified somatic structural variants
|`breakpoints_double.csv` | csv file containing detailed information about identified breakpoints in bam file
|`read_ids.csv`| a csv file containing read ids associated with each identified SV
|`read_qual.txt` | file containing quality statistics about identified segements
|`severus.log` | log file

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
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`germline/sample_GERMLINE_VEP.vcf.gz` | Annotated germline indel and SNV vcf file |
|`germline/sample_GERMLINE_VEP_summary.html` | Visual summary of germline indel and SNV annotations in html format
|`germline/sample_GERMLINE_VEP.vcf.gz.tbi` | Annotated germline indel and SNV vcf index file |
|`somatic/sample_SOMATIC_VEP.vcf.gz` | Annotated somatic indel and SNV vcf file |
|`somatic/sample_SOMATIC_VEP_summary.html` | Visual summary of somatic indel and SNV annotations in html format
|`somatic/sample_SOMATIC_VEP.vcf.gz.tbi` | Annotated somatic indel and SNV vcf index file |
|`SVs/sample_SV_VEP.vcf.gz` | Annotated somatic structural variant vcf file |
|`SVs/sample_SV_VEP_summary.html` | Visual summary of somatic structural variant annotations in html format
|`SVs/sample_SV_VEP.vcf.gz.tbi` | Annotated somatic structural variant vcf index file |

</details>

### wakhan
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



│   ├── solutions_ranks.tsv
```

| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

### pipeline_info

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
| File    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`execution_report_{DATE}.hmtl` | summary of pipeline resource and timing usage in a html report  |
|`execution_timeline_{DATE}.hmtl` | a graphical summary of the timing of each module's task over the course of the pipeline run |
|`lrsomatic_softwar_mqc_versions.yml` | summary of the versions of each tool used by the pipeline|
|`params_{DATE}.json` | summary of the paramaters used in the pipeline|
|`pipeline_dag_{DATE}.html` | flow chart summarizing the pipeline run |

</details>
