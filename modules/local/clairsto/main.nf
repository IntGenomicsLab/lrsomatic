// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process CLAIRSTO {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs-to:v0.3.1':
        'hkubal/clairs-to:v0.3.1' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(tumour_bam)
    path(ref)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_somatic.vcf.gz"), emit: somatic_vcf
    tuple val(meta), path("*_germline.vcf.gz"), emit: germline_vcf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO : change to reflect the meta data information from Laurens
    def platform = task.ext.platform ?: 'ont'
    def output_dir = "${meta.id}_clairs_output"
    def SNV_VCFGZ="${output_dir}/snv.vcf.gz"
    def INDEL_VCFGZ="${output_dir}/indel.vcf.gz"
    def SNV_VCF="${output_dir}/snv.vcf"
    def INDEL_VCF="${output_dir}/indel.vcf"
    def SOMATIC_VCF="${output_dir}/somatic.vcf"
    def NONSOMATIC_VCF="${output_dir}/germline.vcf"

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    /opt/bin/run_clairs_to \\
        --tumor_bam_fn  \\
        --ref_fn ${ref} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir ${output_dir} \\
        --use_heterozygous_snp_in_normal_sample_and_normal_bam_for_intermediate_phasing True \\
        --remove_intermediate_dir \\
        --use_longphase_for_intermediate_phasing True \\
        --use_longphase_for_intermediate_haplotagging True \\
        --conda_prefix /opt/micromamba/envs/clairs-to
   
    # Unzip ;)
    gunzip $SNV_VCFGZ $INDEL_VCFGZ

    # Extract header from snv.vcf (lines starting with '##' or the first line starting with '#')
    awk '/^##/ {print} /^#CHROM/ {print; exit}' "$SNV_VCF" > header.vcf

    # Extract data (excluding headers) from both files
    awk '!/^#/' "$SNV_VCF" > snv_data.tmp
    awk '!/^#/' "$INDEL_VCF" > indel_data.tmp

    # Combine data
    cat snv_data.tmp indel_data.tmp > combined_data.tmp

    # Split into somatic and non-somatic files based on 7th column containing "NonSomatic"
    awk '$7 ~ /NonSomatic/' combined_data.tmp > nonsomatic_data.tmp
    awk '$7 !~ /NonSomatic/' combined_data.tmp > somatic_data.tmp

    # Add header to output files
    cat header.vcf nonsomatic_data.tmp > "$NONSOMATIC_VCF"
    cat header.vcf somatic_data.tmp > "$SOMATIC_VCF"

    # Cleanup temporary files
    rm header.vcf snv_data.tmp indel_data.tmp combined_data.tmp nonsomatic_data.tmp somatic_data.tmp
    
    # Zip the resulting files
    gzip $NONSOMATIC_VCF $SOMATIC_VCF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(/opt/bin/run_clairs --version | sed 's/ClairS version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}_clairs_output"
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/output.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(/opt/bin/run_clairs --version | sed 's/ClairS version: //')
    END_VERSIONS
    """
}
