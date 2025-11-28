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
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

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

### Pipeline options

| Parameter        | Description                                                                                                                                                                  |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-input`         | Full file path to input samplesheet, must be in `.csv` format and conform to specifications noted above                                                                      |
| `--genome`       | Specified genome assembly, support is given for `GRCh38` and `CHM13`                                                                                                         |
| `--normal_fiber` | A boolean which skips fiber-seq processing on normal files (on those which have fiber-seq for the tumor). Default = `true` (_does not skip fiber-seq processing for normal_) |

#### Skipping options:

| Parameter         | Description                                                                                                 |
| ----------------- | ----------------------------------------------------------------------------------------------------------- |
| `--skip_qc`       | A boolean to skip all QC steps, including `mosdepth`, `samtools`,`fibertools`, `cramino`. Default = `false` |
| `--skip_fiber`    | A boolean to skip all `fibertools` related modules. Default = `false`                                       |
| `--skip_cramino`  | A boolean to skip `cramino`. Default = `false`                                                              |
| `--skip_mosdepth` | A boolean to skip `mosdepth`. Default = `false`                                                             |
| `--skip_ascat`    | A boolean to skip `ascat`. Default = `false`                                                                |
| `--skip_bamstats` | A boolean to skip `bamstats`. Default = `false`                                                             |
| `--skip_wakhan`   | A boolean to skip `wakhan`. Default = `false`                                                               |
| `--skip_vep`      | A boolean to skip `vep`. Default = `false`                                                                  |

#### VEP options:

| Parameter             | Description                                                                                                                                      |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--vep_cache`         | Full path to a vep cache. If left blank, this will default to pulling from this [Annotation Cache Storage](https://annotation-cache.github.io/). |
| `--vep_cache_version` | Integer specifying version of vep cache. Default = `113`                                                                                         |
| `--vep_args`          | A string specifying arguments to vep. Default = `"--everything --filter_common --per_gene --total_length --offline --format vcf"`                |
| `--vep_custom`        | A full path to a vcf file containing custom variants for annotation. Must be bgzipped and have `.vcf.gz` format. Default = `null`                |
| `--vep_custom_tbi`    | A full path to a index file for cutom vcf for vep. Default = `null`                                                                              |

#### Minimap2 Options

| Parameter                     | Description                                                                                 |
| ----------------------------- | ------------------------------------------------------------------------------------------- |
| `--minimap2_ont_model`        | specifies which model to use minimap2 with for ONT samples. Default = `null`                |
| `--minimap2_pb_model`         | specifies which model to use minimap2 with for PacBio samples. Default = `null`             |
| `--save_secondary_alignments` | A boolean to specify if secondary alignmetns are kept in aligned bam file. Defualt = `true` |

#### ASCAT Options

| Parameter                     | Description                                                                                       |
| ----------------------------- | ------------------------------------------------------------------------------------------------- |
| `--ascat_ploidy`              | integer to enforce a given ploidy value. Default = `null`                                         |
| `--ascat_purity`              | integer to enforce a given purity value. Default = `null`                                         |
| `--ascat_min_base_qual`       | integer to specify a minimum base quality for ascat's allele counter. Default = `20`              |
| `--ascat_min_counts`          | integer to specify a minimum number of counts for ascat's allele counter. Default = `10`          |
| `--ascat_min_map_qual`        | integer to specify a minimum mapping quality for ascat's allele counter. Default = `10`           |
| `--ascat_penalty`             | integer to specify a penalty value for ascat. Default = `150`                                     |
| `--ascat_longread_bins`       | integer to specify the binsize for ascat long reads. Default = `2000`                             |
| `--ascat_allelecounter_flags` | flags to pass to ascat's allele counter. Default = `"-f 0"`                                       |
| `--ascat_chroms`              | string to enforce a subset of chromosomes on the sample, ie `"(c(1:21,'X','Y')). Default = `null` |

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
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
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

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

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
