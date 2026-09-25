# IntGenomicsLab/lrsomatic: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file of the following form:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber
sample1,tumour.bam,normal.bam,ont,female,n
sample2,tumour.bam,,ont,female,y
sample3,tumour.bam,,pb,male,n
sample4,tumour.bam,normal.bam,pb,male,y
```

lrsomatic extracts information from the bam header files to decide which models to use for Clair3, ClairS, or ClairS-TO. However, this can optionally be specified manually. You can do this for one or many samples, if the field is left blank, the pipeline will default to extracting this information. You can specify this by creating your csv in the following form:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber,clair3_model,clairSTO_model,clairS_model
sample1,tumour.bam,normal.bam,ont,female,n
sample2,tumour.bam,,ont,female,y
sample3,tumour.bam,normal.bam,pb,male,n,r1041_e82_400bps_sup_v420,,ont_r10_dorado_sup_5khz_ssrs
sample4,tumour.bam,normal.bam,pb,male,y
```

Use the `input` parameter to specify the location to this input csv.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber
sample1,tumour1.bam,normal.bam,ont,female,n
sample1,tumour2.bam,,ont,female,y
sample1,tumour3.bam,,pb,male,n
```

### Full Description of Samplesheet Columns

| Column           | Description                                                                                                                                                                            |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`         | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `bam_tumor`      | Full path to BAM file for the tumor. File must end in `.bam`.                                                                                                                          |
| `bam_normal`     | Full path to BAM file for the tumor. File must end in `.bam`.                                                                                                                          |
| `platform`       | A string specifying the platform used for sequencing, can be either `pb` for PacBio sequencing data or `ont` for Oxford Nanopore sequencing data                                       |
| `sex`            | A string specifying the biological sex of the sample, can either be `m` or `f`                                                                                                         |
| `fiber`          | A string specifying if the sample has been subjected to Fiber-seq. Can either be `y` or `n`                                                                                            |
| `clair3_model`   | A string describing which model is to be used for Clair3's small variant calling (_optional_)                                                                                          |
| `clairSTO_model` | A string describing which model is to be used for ClairS-TO's small variant calling (_optional_)                                                                                       |
| `clairS_model`   | A string describing which model is to be used for ClairS's small variant calling (_optional_)                                                                                          |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run IntGenomicsLab/lrsomatic --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.

```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run IntGenomicsLab/lrsomatic -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## CHM13 Support

Our pipeline fully supports CHM13 and most reference and annotation files are automatically downloaded when specifying `--genome CHM13`.

However, VEP will need a bit of additional setup. The VEP cache for CHM13 needs to be manually downloaded. This can be done using the following code. Feel free to change any of the paths, ensuring that the correct path is pointed to in the pipeline parameters.

Download CHM13 Cache:

```bash
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/indexed_vep_cache/Homo_sapiens-GCA_009914755.4-2022_10.tar.gz
tar xzf Homo_sapiens-GCA_009914755.4-2022_10.tar.gz
```

Then you can run the pipeline as follows:

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome CHM13 \
  --vep_cache $HOME/.vep \
  --vep_cache_version 107 \
  -profile docker
```

If you want to run with a CHM13 reference without using `--genome CHM13` (for example, via a custom FASTA or configuration), you must also specify `--vep_genome T2T-CHM13v2.0` and `--vep_species homo_sapiens_gca009914755v4`.

This cache carries gene and transcript models only, so every pathogenicity score on CHM13 has to come from a plugin. With `--genome CHM13` the pipeline fetches those itself; note that only the two predictors keyed on protein rather than genome coordinates can reach the assembly at all. See [VEP plugins](#vep-plugins) for which those are and why.

Setting `--vep_genome T2T-CHM13v2.0` by hand, without `--genome CHM13`, resolves no plugin defaults — supply the `--vep_*` paths yourself, or accept VEP without plugin annotation.

For mutational signatures, `--genome CHM13` selects the `CHM13-T2T` SigProfilerMatrixGenerator genome. Its payload is not on the AlexandrovLab FTP yet, so `--download_sigprofiler_genome` fetches it from the IntGenomicsLab Globus collection (`--sigprofiler_genome_url`); see [Mutational Signature Options](#mutational-signature-options).

For structural variants, the CHM13 panel of normals is a merged panel combining the 1000 Genomes CHM13 panel shipped with SEVERUS and the ASAP cohort, with median confidence intervals per breakpoint. The pipeline exposes it as `--pon_file` and hands it to SEVERUS via that tool's own `--PON` flag; it is downloaded automatically with `--genome CHM13`. GRCh38 continues to use the 1000 Genomes panel shipped with SEVERUS.

For tumour-only small variants, ClairS-TO separates germline from somatic calls with a panel of normals and with its Verdict module, which tags each call as germline, somatic or subclonal somatic from tumour purity and allele-specific copy number. `--genome CHM13` supplies five CHM13 PON VCFs (gnomAD, dbSNP, 1000 Genomes, CoLoRSdb and ASAP), which **replace** the GRCh38 databases inside the container. Unless `--skip_ascat` is set, purity and copy number come from the pipeline's own ASCAT run (`CLAIRSTO_VERDICT_TAG`); only with `--skip_ascat` does ClairS-TO estimate them itself, from assembly-specific loci, allele and GC content files. A GRCh38 resource set on a CHM13 run leaves germline variants untagged.

With `--genome CHM13 --skip_ascat` the pipeline builds a CHM13 resource set from the ASCAT files it already downloads, so no extra setup is needed. LogR correction is GC-only, as ClairS-TO recommends for CHM13: no replication timing file is published for the assembly. Without `--skip_ascat` nothing is built, because the tagging comes from ASCAT's own tables.

To use a resource set of your own — another assembly, or a CHM13 set carrying an `RT_<name>.txt` for replication timing correction — pass `--clairsto_cna_resources` together with `--skip_ascat`. Without `--skip_ascat` it is ignored, with a warning. Expected layout:

```
<dir>/
  loci_files/<prefix>chr1.txt ... <prefix>chr22.txt, <prefix>chrX.txt
  allele_files/<prefix>chr1.txt ... <prefix>chr22.txt, <prefix>chrX.txt
  GC_<name>.txt
  RT_<name>.txt        (optional; without it, correction is GC-only)
```

Keep one resource set per directory: ClairS-TO derives the per-contig prefix from the single `*chr1.txt` and takes exactly one `GC_*.txt`. The directory must be self-contained, since only the directory itself is staged into the container; materialise symlinked content with `cp -rL`. Both are checked at startup.

If the loci cannot belong to the reference, ClairS-TO disables Verdict with a warning rather than applying the wrong coordinates. On a `--skip_ascat` run, look for `VERDICT CNA RESOURCE DIRECTORY` in the ClairS-TO log to confirm which set was used; otherwise Verdict runs outside ClairS-TO and the tables the tags came from are published next to the VCFs.

### Pipeline options

| Parameter        | Description                                                                                                                                                                  |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-input`         | Full file path to input samplesheet, must be in `.csv` format and conform to specifications noted above                                                                      |
| `--genome`       | Specified genome assembly, support is given for `GRCh38` and `CHM13`                                                                                                         |
| `--normal_fiber` | A boolean which skips fiber-seq processing on normal files (on those which have fiber-seq for the tumor). Default = `true` (_does not skip fiber-seq processing for normal_) |

#### Skipping options:

| Parameter              | Description                                                                                                                                                                                                                                                             |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--skip_qc`            | A boolean to skip all QC steps, including `mosdepth`, `samtools`,`fibertools`, `cramino`. Default = `false`                                                                                                                                                             |
| `--skip_fiber`         | A boolean to skip all `fibertools` related modules. Default = `false`                                                                                                                                                                                                   |
| `--skip_cramino`       | A boolean to skip `cramino`. Default = `false`                                                                                                                                                                                                                          |
| `--skip_mosdepth`      | A boolean to skip `mosdepth`. Default = `false`                                                                                                                                                                                                                         |
| `--skip_ascat`         | A boolean to skip `ascat`. ClairS-TO's Verdict germline tagging then falls back to Verdict's own purity and copy number estimate, which is still up to 0.14 from ASCAT's on the samples it was measured on — see [Verdict tags](output.md#clairs-to). Default = `false` |
| `--skip_savana`        | A boolean to skip `savana` (SV + copy-number calling). Default = `false`                                                                                                                                                                                                |
| `--skip_bamstats`      | A boolean to skip `bamstats`. Default = `false`                                                                                                                                                                                                                         |
| `--skip_wakhan`        | A boolean to skip `wakhan`. Default = `false`                                                                                                                                                                                                                           |
| `--skip_vep`           | A boolean to skip `vep`. Default = `false`                                                                                                                                                                                                                              |
| `--skip_m6a`           | A boolean to skip `fibertools_m6a`, used if you have m6a calls but would still like nucleosome positions for PacBio data (ONT data is required to have m6a calls). Default = `false`                                                                                    |
| `--skip_nanoplot`      | A boolean to skip NanoPlot QC on aligned and unaligned BAM files. Default = `false`                                                                                                                                                                                     |
| `--skip_normalfiber`   | A boolean to skip fibertools processing for the normal sample. Default = `false`                                                                                                                                                                                        |
| `--skip_modcall`       | A boolean to skip Longphase `modcall`, the 5mC base-modification calling whose VCF is used as extra evidence during phasing. Unrelated to the modkit pileup (see `--skip_modkit`). Default = `false`                                                                    |
| `--skip_modkit`        | A boolean to skip the modkit pileup step. Default = `false`                                                                                                                                                                                                             |
| `--skip_whatshapstats` | A boolean to skip WhatsHap phasing statistics. Default = `false`                                                                                                                                                                                                        |
| `--skip_signatures`    | A boolean to skip mutational signature analysis (SigProfilerMatrixGenerator + SigProfilerAssignment). Default = `false`                                                                                                                                                 |
| `--skip_report`        | A boolean to skip the final per-sample HTML report. Default = `false`                                                                                                                                                                                                   |

#### Modkit options:

| Parameter         | Description                                                                                                                                                                                                                                                                                                                        |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--modkit_args`   | Arguments passed to `modkit pileup`. The value replaces the default rather than extending it, so repeat `--cpg --modified-bases 5mC` when adding flags. An empty value gives the unfiltered pileup: use the `=` form, `--modkit_args=''`, or an empty `modkit_args` entry in a params file. Default = `--cpg --modified-bases 5mC` |
| `--modkit_phased` | Run `modkit pileup --phased` on the Longphase-haplotagged BAMs, producing `_hp1`, `_hp2` and `_combined` bedMethyl files per sample. The pileup then requires small-variant calling, phasing and haplotagging to have succeeded for the sample. Default = `false`                                                                  |

The pileup runs a patched modkit 0.6.4 image (`ghcr.io/ljwharbers/modkit`, `linux/amd64` only, built from [ljwharbers/modkit@pacbio-conflict-fix](https://github.com/ljwharbers/modkit/tree/pacbio-conflict-fix)): stock modkit drops PacBio HiFi reads whose 5mC and 5hmC probabilities sum above 1 and returns empty `--cpg` pileups on them ([nanoporetech/modkit#612](https://github.com/nanoporetech/modkit/issues/612)). Conda and `arm64` are not supported for this step; use `--skip_modkit` there.

#### LONGPHASE options:

| Parameter                       | Description                                                                        |
| ------------------------------- | ---------------------------------------------------------------------------------- |
| `--longphase_tag_supplementary` | Include supplementary alignments in Longphase haplotype tagging. Default = `false` |

#### VEP options:

| Parameter              | Description                                                                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--vep_cache`          | Full path to a vep cache. If left blank, this will default to pulling from this [Annotation Cache Storage](https://annotation-cache.github.io/). |
| `--vep_cache_version`  | Integer specifying version of vep cache. Default = `113`                                                                                         |
| `--vep_args`           | A string specifying arguments to vep. Default = `"--everything --filter_common --per_gene --total_length --offline --format vcf --vcf"`          |
| `--vep_custom`         | A full path to a vcf file containing custom variants for annotation. Must be bgzipped and have `.vcf.gz` format. Default = `null`                |
| `--vep_custom_tbi`     | A full path to a index file for cutom vcf for vep. Default = `null`                                                                              |
| `--download_vep_cache` | A boolean to automatically download the VEP cache if not found locally. Default = `false`                                                        |

#### VEP plugin options:

The plugin data is **on by default** on `--genome GRCh38` and `--genome CHM13`; setting a
parameter overrides its default, `--skip_vep_plugins` turns the set off, and CADD and EVE are
opt-in. See [VEP plugins](#vep-plugins) for sizes, licence terms and per-assembly availability.

| Parameter                    | Description                                                                                                              |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `--skip_vep_plugins`         | Annotate with VEP alone, fetching no plugin data. Default = `false`                                                      |
| `--vep_alphamissense`        | AlphaMissense GRCh38 score file, for the `AlphaMissense` plugin. GRCh38 only                                             |
| `--vep_alphamissense_tbi`    | Index for `--vep_alphamissense`. Required whenever `--vep_alphamissense` is set                                          |
| `--vep_alphamissense_aa`     | AlphaMissense protein-space release, or a table already built from it, for the `AlphaMissenseProtein` plugin. CHM13 only |
| `--vep_alphamissense_aa_tbi` | Index for `--vep_alphamissense_aa`. Required whenever `--vep_alphamissense_aa` is set                                    |
| `--vep_polyphen_sift_db`     | Ensembl pangenome PolyPhen/SIFT SQLite database, for the `PolyPhen_SIFT` plugin. Needed on CHM13 only                    |
| `--vep_clinvar`              | ClinVar VCF, added as a VEP `--custom` annotation                                                                        |
| `--vep_clinvar_tbi`          | Index for `--vep_clinvar`. Required whenever `--vep_clinvar` is set                                                      |
| `--vep_clinvar_md5`          | Expected MD5 of a remote `--vep_clinvar`, checked after download. Dropped when `--vep_clinvar` is overridden             |
| `--vep_clinvar_tbi_md5`      | Expected MD5 of the downloaded `--vep_clinvar_tbi`. Set on CHM13 only; dropped when either ClinVar file is overridden    |
| `--vep_clinvar_fields`       | Comma-separated ClinVar INFO fields to carry through. Default = `"CLNSIG,CLNREVSTAT,CLNDN"`                              |
| `--vep_cadd_snv`             | CADD SNV score file, for the `CADD` plugin. No default — 81 GB, so opt-in; prefer a local path. GRCh38 only              |
| `--vep_cadd_snv_tbi`         | Index for `--vep_cadd_snv`. Required whenever `--vep_cadd_snv` is set                                                    |
| `--vep_cadd_indel`           | CADD indel score file, for the `CADD` plugin. No default — opt-in. GRCh38 only                                           |
| `--vep_cadd_indel_tbi`       | Index for `--vep_cadd_indel`. Required whenever `--vep_cadd_indel` is set                                                |
| `--vep_revel`                | REVEL release zip, or a prepared score file, for the `REVEL` plugin. GRCh38 only                                         |
| `--vep_revel_tbi`            | Index for `--vep_revel`. Required only when `--vep_revel` is an already-prepared file                                    |
| `--vep_eve`                  | EVE release zip, or a merged VCF, for the `EVE` plugin. **Not** on by default. GRCh38 only. Default = `null`             |
| `--vep_eve_tbi`              | Index for `--vep_eve`. Required only when `--vep_eve` is an already-merged VCF. Default = `null`                         |

#### Minimap2 Options

| Parameter                    | Description                                                                                 |
| ---------------------------- | ------------------------------------------------------------------------------------------- |
| `--minimap2_ont_model`       | specifies which model to use minimap2 with for ONT samples. Default = `null`                |
| `--minimap2_pb_model`        | specifies which model to use minimap2 with for PacBio samples. Default = `null`             |
| `--save_secondary_alignment` | A boolean to specify if secondary alignments are kept in aligned bam file. Default = `true` |

#### ASCAT Options

| Parameter                     | Description                                                                                                                                                                                                                 |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--ascat_ploidy`              | integer to enforce a given ploidy value. Default = `null`                                                                                                                                                                   |
| `--ascat_purity`              | integer to enforce a given purity value. Default = `null`                                                                                                                                                                   |
| `--ascat_min_base_qual`       | integer to specify a minimum base quality for ascat's allele counter. Default = `20`                                                                                                                                        |
| `--ascat_min_counts`          | integer to specify a minimum number of counts for ascat's allele counter. Default = `10`                                                                                                                                    |
| `--ascat_min_map_qual`        | integer to specify a minimum mapping quality for ascat's allele counter. Default = `10`                                                                                                                                     |
| `--ascat_penalty`             | integer to specify a penalty value for ascat. Default = `150`                                                                                                                                                               |
| `--ascat_longread_bins`       | integer to specify the binsize for ascat long reads. Default = `2000`                                                                                                                                                       |
| `--ascat_allelecounter_flags` | flags to pass to ascat's allele counter. Default = `"-f 0"`                                                                                                                                                                 |
| `--ascat_chroms`              | string to enforce a subset of chromosomes on the sample, ie `"(c(1:21,'X','Y')). Default = `null`                                                                                                                           |
| `--ascat_allele_files`        | A full path to a zipped folder containing allele files for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Must be zipped and have `.zip` format. Default = `null`                             |
| `--ascat_loci_files`          | A full path to a zipped folder containing loci files for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Must be zipped and have `.zip` format. Default = `null`                               |
| `--ascat_gc_file`             | A full path to a GC correction file for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Optionally can be zipped and have either `.txt` or `.txt.zip` format. Default = `null`                 |
| `--ascat_rt_file`             | A full path to a replication timing correction file for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Optionally can be zipped and have either `.txt` or `.txt.zip` format. Default = `null` |
| `--ascat_pdf_plots`           | string to enable output pltos in pdf format. Default = `false`                                                                                                                                                              |

#### Fibertools Options

| Parameter           | Description                                                                      |
| ------------------- | -------------------------------------------------------------------------------- |
| `--autocorrelation` | A boolean to enable autocorrelation computation in fibertools. Default = `false` |

#### SEVERUS Options

| Parameter              | Description                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------ |
| `--severus_minsupport` | Minimum number of supporting reads required for SEVERUS to call an SV. Default = `3` |

#### SAVANA Options

SAVANA reuses the haplotagged BAMs from `PHASING_HAPLOTYPING` and the same reference FASTA/index as
every other caller. Matched tumor/normal samples run `savana run` + `savana classify` + `savana cna`,
using the phased germline VCF as the SNP allele-frequency source for copy-number fitting; tumor-only
samples run `savana to` instead (which chains the equivalent steps internally), using SAVANA's bundled
1000 Genomes population SNP set rather than the tumour's own calls. Both modes are restricted to
canonical chromosomes via a genome-specific `--contigs` file.

| Parameter                | Description                                                                                                           |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `--savana_pb_minsupport` | Minimum supporting reads for SAVANA to call a variant on PacBio samples (`--min_support` with `--pb`). Default = `10` |

#### Report Options

| Parameter             | Description                                                                                                                                                                                                                                  |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--report_gene_panel` | Gene panel(s) applied when the report opens, as a comma-separated list. Each entry is `none` (no filtering), a builtin panel name (`lymphoid` or `sarcoma`), or a path to a TSV file with a `gene` column. Default = `null`, i.e. unfiltered |

The report is rendered by [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) from
`ghcr.io/ljwharbers/lrsomatic-report` (`ghcr.io/ljwharbers/lrsomatic-report-sif` under
Singularity/Apptainer), pinned to a tag in
[`modules/local/lrsomaticreport/main.nf`](../modules/local/lrsomaticreport/main.nf). Images are
built for `linux/amd64` only and **Conda is not supported for this step**; use `--skip_report`
under `-profile conda`/`mamba`. Builtin gene panels live in
[`assets/gene_lists/`](../assets/gene_lists/README.md): a TSV dropped there becomes a valid
`--report_gene_panel` value.

Gene panel filtering is a view, not a filter on the data: every builtin panel is embedded in the
report and can be ticked or unticked in the browser; `--report_gene_panel` only decides which are
ticked on load. A custom panel is a tab-separated file with a header row and at least a `gene`
column:

```tsv
gene	panel	note
TP53	mypanel	Tumour suppressor
KRAS	mypanel	Oncogene
```

Optional `chrom`, `start` and `end` columns (all three or none, with the reference declared as a
`# reference: hg38` line or a `reference` column) switch SV matching from the VEP gene symbol to
position, and an optional `applies_to` column (`snv`, `sv` or blank) scopes a gene to one table.
See [`assets/gene_lists/README.md`](../assets/gene_lists/README.md) for the full format.

> **The builtin `lymphoid` panel was rebuilt** in the release that added `applies_to`: 234 rows
> scoped per table, where it had 72 genes filtering both. `--report_gene_panel lymphoid` therefore
> gives different filtered tables across that change. `sarcoma` is unchanged.

```bash
nextflow run IntGenomicsLab/lrsomatic \
    -profile <docker/singularity> \
    --input samplesheet.csv \
    --outdir results \
    --report_gene_panel /path/to/mypanel.tsv
```

##### Applying several panels at once

Pass a comma-separated list to open the report with several panels applied. They are
**unioned**: a variant or SV is kept if it hits any of them. Builtin names and custom paths
can be mixed freely.

```bash
    --report_gene_panel lymphoid,sarcoma
    --report_gene_panel 'lymphoid,/path/to/mypanel.tsv'
```

With two or more panels active, each `panel_hit` entry in the SV table gains a trailing `[panel]`
naming which one matched.

An entry containing a `/` or ending in `.tsv` is read as a panel file; anything else is a builtin
name, so `--report_gene_panel mypanel` is looked up as a builtin even if a file called `mypanel`
exists.

Checked at launch:

- `none` cannot be combined with a real panel.
- A panel file that does not exist, or a builtin name that is not bundled, is an error.
- Two panel files sharing a base name cannot be combined; rename one.

#### WAKHAN Options

| Parameter         | Description                                                                                             |
| ----------------- | ------------------------------------------------------------------------------------------------------- |
| `--wakhan_chroms` | A string specifying a subset of chromosomes for WAKHAN to process, e.g. `"chr1,chr2"`. Default = `null` |

#### Mutational Signature Options

Mutational signature analysis runs [SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) on the PASS SNVs and indels of each sample's phased somatic VCF and fits COSMIC signatures per sample with [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment). It needs SigProfilerMatrixGenerator's per-genome payload (~3 GB), which is not shipped with the pipeline. Either:

- run once with `--download_sigprofiler_genome`; the payload is installed, checksum-verified and published to `<outdir>/cache/sigprofiler/volume`, or
- point `--sigprofiler_genome_dir` at an existing SigProfilerMatrixGenerator volume (a directory containing `tsb/<genome>/`, e.g. one created with `SigProfilerMatrixGenerator install GRCh38 --volume <dir>` or the published cache from a previous run).

Running with neither, and without `--skip_signatures`, stops the pipeline at start-up.

| Parameter                                   | Description                                                                                                                                                                                         |
| ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--sigprofiler_genome_dir`                  | Full path to a SigProfilerMatrixGenerator volume containing `tsb/<sigprofiler_genome>/`. Default = `null`                                                                                           |
| `--download_sigprofiler_genome`             | A boolean to download and install the genome payload during the run (published to `<outdir>/cache/sigprofiler/volume`). Default = `false`                                                           |
| `--sigprofiler_cosmic_version`              | COSMIC reference signature version fitted by SigProfilerAssignment. Default = `3.6`                                                                                                                 |
| `--sigprofiler_exclude_signature_subgroups` | Comma-separated SigProfilerAssignment signature subgroups to exclude from the fit, e.g. `"Artifact_signatures,Lymphoid_signatures"` (see the SigProfilerAssignment documentation). Default = `null` |
| `--sigprofiler_matrix_args`                 | Extra arguments for `SigProfilerMatrixGenerator matrix_generator`. Default = `"--plot"`                                                                                                             |
| `--sigprofiler_assignment_args`             | Extra arguments for `SigProfilerAssignment cosmic_fit`, e.g. `"--make_plots False"`. Default = `null`                                                                                               |

Both tools run from `ghcr.io/ljwharbers/sigprofiler`, which adds CHM13 support not yet in a SigProfiler release: SigProfilerMatrixGenerator from the branch behind [SigProfilerSuite/SigProfilerMatrixGenerator#250](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/pull/250) and SigProfilerAssignment from [ljwharbers/SigProfilerAssignment](https://github.com/ljwharbers/SigProfilerAssignment/tree/chm13-t2t-support), with COSMIC SBS/DBS signatures renormalised to CHM13. ID83 signatures always use the GRCh37 set. Conda is not supported for this step.

#### Variant Filtering and Combining Options

These options control how variants from multiple callers are filtered and merged.

| Parameter                      | Description                                                                                         |
| ------------------------------ | --------------------------------------------------------------------------------------------------- |
| `--germline_var_keep`          | Expression or threshold for retaining germline variants after calling. Default = `null`             |
| `--somatic_var_keep`           | Expression or threshold for retaining somatic variants after calling. Default = `null`              |
| `--germline_var_combine`       | Strategy for combining germline variant caller outputs (e.g. union, intersection). Default = `null` |
| `--somatic_var_combine`        | Strategy for combining somatic variant caller outputs (e.g. union, intersection). Default = `null`  |
| `--prioritize_caller_germline` | Comma-separated caller priority order used when combining germline calls. Default = `null`          |
| `--prioritize_caller_somatic`  | Comma-separated caller priority order used when combining somatic calls. Default = `null`           |

#### PON Options

| Parameter                | Description                                                                                                                                                                                                                                                 |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--clairsto_pon_vcfs`    | Full path to one or more Panel of Normals VCF files for ClairS-TO small variant filtering. Default = `null`, meaning the set for `--genome` is used: four VCFs on GRCh38, five on CHM13. Supplying VCFs replaces that set entirely rather than adding to it |
| `--clairsto_pon_flags`   | Population allele matching flags for ClairS-TO PON VCFs (one per VCF, comma-separated). Default = `null`, meaning the flags that go with the `--genome` set. There must be exactly one flag per VCF, or the run stops at startup                            |
| `--deepsomatic_pon_vcfs` | Full path to one or more bgzipped, tabix-indexed PON VCF files (for example, `.vcf.gz`) passed to DeepSomatic `--population_vcfs`. If not set, uses container-bundled defaults in tumor-only mode or no PON in paired mode. Default = `null`                |

#### Advanced Options

| Parameter         | Description                                                                                                                                                                 |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--use_gpu`       | A boolean to enable GPU acceleration for DeepVariant and DeepSomatic. Requires a GPU-enabled compute environment. Default = `false`                                         |
| `--generate_gvcf` | A boolean to enable gVCF output from DeepVariant (germline) and DeepSomatic (somatic). gVCF files include calls at all positions, not just variant sites. Default = `false` |

#### Genome-Derived Parameters

The following parameters are automatically populated from the `--genome` iGenomes configuration and do not normally need to be set manually. They can be overridden when using a custom genome or reference build not present in the iGenomes configuration.

| Parameter                                                                          | Description                                                                                                                                                                                                                                         |
| ---------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--fasta`                                                                          | Full path to the reference FASTA file. Auto-populated from `--genome`. Override for custom genomes.                                                                                                                                                 |
| `--bed_file`                                                                       | BED file of callable/target regions passed to SEVERUS for SV calling. Auto-populated from `--genome`.                                                                                                                                               |
| `--pon_file`                                                                       | Panel of Normals breakpoint table (bgzipped CSV) for SEVERUS somatic SV filtering in tumor-only mode. Auto-populated from `--genome`.                                                                                                               |
| `--centromere_bed`                                                                 | BED file of centromere coordinates passed to WAKHAN. Auto-populated from `--genome`.                                                                                                                                                                |
| `--genome_name`                                                                    | Assembly name string passed to ASCAT for genome-specific reference file selection. Auto-populated from `--genome`.                                                                                                                                  |
| `--vep_genome`                                                                     | VEP genome identifier (e.g. `GRCh38`, `T2T-CHM13v2.0`). Auto-populated from `--genome`. Override for CHM13 or custom assemblies.                                                                                                                    |
| `--vep_species`                                                                    | VEP species identifier. Auto-populated from `--genome`. Override for non-standard assemblies (e.g. `homo_sapiens_gca009914755v4` for CHM13).                                                                                                        |
| `--sigprofiler_genome`                                                             | SigProfilerMatrixGenerator / SigProfilerAssignment genome name (`GRCh38` or `CHM13-T2T`). Auto-populated from `--genome`.                                                                                                                           |
| `--sigprofiler_genome_url`                                                         | URL of the `<sigprofiler_genome>.tar.gz` payload for `--download_sigprofiler_genome`; unset means the AlexandrovLab FTP. Auto-populated from `--genome` (set for CHM13 only).                                                                       |
| `--ascat_allele_files`, `--ascat_loci_files`, `--ascat_gc_file`, `--ascat_rt_file` | ASCAT loci, allele, GC content and replication timing files. Auto-populated from `--genome`. Also used to build the ClairS-TO Verdict CNA resources on CHM13; CHM13 has no RT file, so that correction is GC-only.                                  |
| ClairS-TO PON VCFs                                                                 | The four (GRCh38) or five (CHM13) VCFs behind `--clairsto_pon_vcfs`: `gnomad`, `dbsnp`, `onekgenomes`, `colors`, and on CHM13 also `asap`. Auto-populated from `--genome`; override with `--clairsto_pon_vcfs` and `--clairsto_pon_flags` together. |

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull IntGenomicsLab/lrsomatic
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [IntGenomicsLab/lrsomatic releases page](https://github.com/IntGenomicsLab/lrsomatic/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## VEP plugins

VEP is given extra pathogenicity and clinical-significance annotation through plugins. On
`--genome GRCh38` and `--genome CHM13` these are **on by default**: the pipeline downloads the
resources for the assembly and reshapes the ones VEP cannot read as published. Pass
`--skip_vep_plugins` to annotate with VEP alone; the plugins are also off under `--igenomes_ignore`
or with any other `--genome`.

Every default points at the resource's original source, except the AlphaMissense GRCh38 tabix index
and the CHM13 protein-space table, which the lab hosts because neither is published in a form VEP
can use (both CC BY 4.0 — see `CITATIONS.md`). Complying with each licence remains your
responsibility, and **CADD, REVEL and EVE are free for non-commercial use only** — see the callout
at the end of this section.

Plugins are applied to the germline and somatic VEP runs, not to the structural-variant run.

### What is enabled by default, per assembly

| Tool              | GRCh38                                   | CHM13                         | Default download                                    |
| ----------------- | ---------------------------------------- | ----------------------------- | --------------------------------------------------- |
| **SIFT**          | already in the VEP cache, no file needed | `PolyPhen_SIFT` plugin        | CHM13: 13 GB database                               |
| **PolyPhen**      | already in the VEP cache, no file needed | `PolyPhen_SIFT` plugin        | as above, the same database                         |
| **AlphaMissense** | `AlphaMissense` plugin                   | `AlphaMissenseProtein` plugin | 613 MB (GRCh38) / 1.1 GB (CHM13)                    |
| **ClinVar**       | `--custom` annotation                    | `--custom`, CHM13-lifted VCF  | 193 MB (GRCh38) / 99 MB (CHM13)                     |
| **CADD**          | `CADD` plugin, opt-in                    | not available — see below     | none — `--vep_cadd_snv` enables it (81 GB + 1.2 GB) |
| **REVEL**         | `REVEL` plugin                           | not available — see below     | 667 MB release zip                                  |
| **EVE**           | `EVE` plugin, opt-in                     | not available — see below     | none — `--vep_eve` enables it (9.6 GB)              |

That comes to roughly **1.4 GB on GRCh38** and **14.3 GB on CHM13** for a first run, the CHM13
figure being almost entirely the pangenome database. On GRCh38, SIFT and PolyPhen already come out
of the VEP cache through the default `--everything`, so `--vep_polyphen_sift_db` is only used on
CHM13.

Two resources are left opt-in, both because of their size:

- **CADD** — the SNV table alone is 81 GB, plus 1.2 GB of indels, and it is re-staged on every run
  that is not a `-resume`. Enable it with `--vep_cadd_snv <file> --vep_cadd_snv_tbi <file>` and, for
  indels, `--vep_cadd_indel` with its `_tbi`; download once and point at a local copy.
- **EVE** — a 9.6 GB zip of per-protein VCFs that then has to be merged, for a predictor
  AlphaMissense largely covers. Enable it with `--vep_eve https://evemodel.org/api/proteins/bulk/download/`
  (or a path to the zip, or to a VCF you have already merged).

### Downloading and reusing the data

Prepared REVEL and EVE files are published to `<outdir>/vep_plugins/`, alongside their indexes, so
a later run can point `--vep_revel` / `--vep_revel_tbi` (or the `--vep_eve` pair) at them and skip
both the download and the reshaping. The two AlphaMissense tables are fetched already indexed.

> [!IMPORTANT]
> `<outdir>/vep_plugins/` holds data that is free for non-commercial use only. Exclude it when you
> share or archive a results directory — passing it on is redistribution, which those licences do
> not grant you.

Within a run, downloaded resources are shared by every VEP task. Across runs they are reused only
with `-resume`; a fresh run in the same work directory downloads everything again.

For repeat runs, download once and pass local paths, which skips both the download and the
preparation step:

```bash
# once
bin/prepare_vep_plugin_data.sh revel /data/vep/revel-v1.3_all_chromosomes.zip /data/vep

# then, on every run
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv --outdir ./results --genome GRCh38 \
  --vep_cadd_snv /data/vep/whole_genome_SNVs.tsv.gz \
  --vep_cadd_snv_tbi /data/vep/whole_genome_SNVs.tsv.gz.tbi \
  --vep_revel /data/vep/revel_grch38.tsv.gz \
  --vep_revel_tbi /data/vep/revel_grch38.tsv.gz.tbi \
  -profile singularity
```

Everything except REVEL and EVE needs no preparation — download it with `curl` and pass the path.
`bin/prepare_vep_plugin_data.sh` is what the pipeline itself runs for those two, so a file you build
with it is identical to one it would prepare; run it with `--help` for usage.

Whether an index is required depends on the shape of what you supply:

- **AlphaMissense, ClinVar and CADD** are used exactly as given, so their `_tbi` parameter is
  always required alongside them. Overriding a data file drops the default index: supply both, or
  neither.
- **A remote ClinVar** (http, https or ftp) is downloaded once per run by `VEPPLUGIN_CLINVAR`, with
  its index, rather than staged by `GERMLINE_VEP` and `SOMATIC_VEP` for every sample: NCBI answers
  the burst of requests a multi-sample run sends with HTTP 503. The index must then be a URL too. The
  download is checked against `--vep_clinvar_md5` (and `--vep_clinvar_tbi_md5` where set), so a
  release re-published under the same name fails the run instead of changing the annotation. Both
  defaults carry the MD5 their host publishes, and the CHM13 default also pins its index; NCBI
  publishes no index checksum. With your own URL, pass its MD5 too, or the pipeline warns that the
  release is not verified. Local and cloud-storage (`s3://`, `gs://`, `az://`) paths are staged as
  given. The files are published to `<outdir>/vep_plugins/` for reuse as local paths.
- **REVEL and EVE** ship as zip archives. Pass a `.zip` and the pipeline unpacks and reshapes it;
  pass a prepared file and its index to use it directly. Remote zips are fetched with `wget` rather
  than staged by Nextflow, and EVE's 9.6 GB archive serves slowly, so expect hours.

### Where each default comes from

| Resource                  | Source                                                                                                                            | Size   | Prepared by the pipeline                   | Licence                      |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------ | ------------------------------------------ | ---------------------------- |
| AlphaMissense (GRCh38)    | `https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz`, with an index we host                                | 613 MB | no, the index is fetched ready             | CC BY 4.0                    |
| AlphaMissense (protein)   | a gene-symbol-keyed table we host, built from the AlphaMissense protein-space release                                             | 1.1 GB | no, fetched ready — see `CITATIONS.md`     | CC BY 4.0                    |
| Pangenome PolyPhen/SIFT   | `https://ftp.ensembl.org/pub/release-115/variation/pangenomes/Human/homo_sapiens_pangenome_PolyPhen_SIFT_20240502.db`             | 13 GB  | no, it is an SQLite database               | Ensembl / EMBL-EBI open      |
| ClinVar (GRCh38)          | `clinvar_20260829.vcf.gz` under `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/` (+ `.tbi`)                | 193 MB | downloaded once, MD5-checked               | public domain                |
| ClinVar (CHM13)           | `clinvar_20240624_GCA_009914755.4.vcf.gz` under `https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/` | 99 MB  | downloaded once, MD5-checked               | public domain                |
| CADD v1.7 SNVs (opt-in)   | `https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz` (+ `.tbi`)                                 | 81 GB  | no, `.tbi` is published                    | free for non-commercial use  |
| CADD v1.7 indels (opt-in) | `https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz` (+ `.tbi`)                         | 1.2 GB | no, `.tbi` is published                    | free for non-commercial use  |
| REVEL v1.3                | `https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip`                                                                | 667 MB | unpacked, re-sorted on GRCh38, indexed     | free for non-commercial use  |
| EVE (opt-in)              | `https://evemodel.org/api/proteins/bulk/download/`                                                                                | 9.6 GB | unpacked, per-protein VCFs merged, indexed | non-commercial; see the site |

### Example: GRCh38

Defaults do the work, so nothing plugin-related is needed:

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  -profile docker
```

### Example: CHM13

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome CHM13 \
  --vep_cache $HOME/.vep \
  --vep_cache_version 107 \
  -profile docker
```

Passing a GRCh38-only resource together with a CHM13 cache is rejected before the run starts,
rather than failing hours later inside VEP or after a long download.

### How annotation reaches CHM13 at all

The CHM13 cache is the Ensembl rapid-release HPRC cache: gene and transcript models only, with no
variation, RefSeq or regulatory data. Every genome-coordinate-keyed score resource is published for
GRCh37/GRCh38 alone, and there are no lifted score tracks for the assembly elsewhere either.

Two of these predictors get there anyway, because they score _proteins_ rather than genome positions:

- **`PolyPhen_SIFT`** keys its lookup on the MD5 of the peptide sequence, which is why Ensembl's
  pangenome database covers the HPRC assemblies.
- **`AlphaMissenseProtein`** (in `assets/vep_plugins/`) keys on gene symbol plus amino-acid
  substitution, using a table built from AlphaMissense's protein-space release. A row is used only
  when both amino acids match what VEP computed for the CHM13 transcript; otherwise it reports
  `aa_mismatch` and no score. The `AlphaMissenseProtein_match` values are listed under
  [plugin fields in `CSQ`](output.md#plugin-fields-in-the-csq-annotation).

### What is not available, and why

- **MutationTaster** — no VEP plugin exists, and the dbNSFP route is academic-only.
- **CADD on CHM13** — CADD scores non-coding positions too, so it has no protein-space form to fall
  back on.
- **REVEL and EVE on CHM13** — neither licence clearly permits distributing a derived table;
  AlphaMissense covers the same class of variant and is CC BY 4.0.
- **SpliceAI** — not currently wired up on either assembly.

> [!IMPORTANT]
> **REVEL is enabled by default and is free for non-commercial use only**; CADD and EVE, if you
> enable them, are the same. The pipeline cannot accept those terms on your behalf: if your work is
> commercial, pass `--skip_vep_plugins`, or set only the resources you are licensed for.

> [!IMPORTANT]
> Check that contig naming matches your reference: the GRCh38 reference used here is GATK-style
> (`chr1`), while NCBI's ClinVar VCF ships Ensembl-style names (`1`). VEP usually reconciles the two
> for `--custom` files, but a name it cannot map annotates nothing rather than raising an error, so
> confirm that `CLNSIG` values appear in an annotated VCF before trusting them.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Large container images under Singularity/Apptainer

The ClairS-TO image is the pipeline's largest custom image (about 3.3 GB as a SIF). Under
Singularity and Apptainer it is pulled as a prebuilt SIF from Docker Hub,
`oras://docker.io/ljwharbers/clairs-to-sif:<tag>`; every other engine runs
`docker.io/ljwharbers/clairs-to:<tag>`. It lives on Docker Hub rather than on ghcr with the
pipeline's other custom images because ghcr redirects each download to a signed URL that expires
at the next 5-minute mark and cuts a stream that is still open then, and Apptainer cannot resume a
cut download: the SIF failed with `PROTOCOL_ERROR` on any link slower than about 10 MB/s. Docker
Hub's download URLs are valid for 50 minutes and are only checked when the request starts, so a
slow pull completes.

What this means in practice:

- On a slow link the pull can exceed a site's `pullTimeout` (Nextflow's default is 20 minutes).
  Raise it with `-c`:

  ```groovy
  singularity.pullTimeout = '2h' // apptainer.pullTimeout under -profile apptainer
  ```

- Docker Hub limits anonymous pulls per source IP address, which a cluster's shared outbound
  address can exhaust. If pulls fail with `TOOMANYREQUESTS`, authenticate with a free Docker Hub
  account by setting `APPTAINER_DOCKER_USERNAME` and `APPTAINER_DOCKER_PASSWORD`
  (`SINGULARITY_DOCKER_USERNAME` / `SINGULARITY_DOCKER_PASSWORD` under Singularity) in the
  launching shell, or with `apptainer remote login --username <user> oras://docker.io`.
- Pulling once into a shared cache saves every later run the download. The file name is the one
  Nextflow derives from the image name:

  ```bash
  apptainer pull "$NXF_SINGULARITY_CACHEDIR/docker.io-ljwharbers-clairs-to-sif-0.5.1-verdict-chm13-c0687e8-flat.img" \
      oras://docker.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13-c0687e8-flat
  ```

  `NXF_SINGULARITY_LIBRARYDIR` is checked first and never written to, so a read-only shared
  directory of prebuilt images works as well.

The other custom SIFs (modkit, SigProfiler and the report image, the largest at 0.8 GB) are still
`oras://ghcr.io/ljwharbers/<name>-sif:<tag>` and can hit the same cut on a very slow link. If one
does, the SIF bytes are the registry blob, which ghcr serves with Range support. Download it
resumably into the cache under the name Nextflow expects (`ghcr.io-ljwharbers-<name>-sif-<tag>.img`)
with `curl -C -`, re-requesting the anonymous token from
`https://ghcr.io/token?scope=repository:ljwharbers/<name>-sif:pull` on each attempt, until the file
reaches the layer size in the image manifest; then check its `sha256sum` against the layer digest
and `chmod +x` it.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
