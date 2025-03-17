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

process CLAIRS {
    tag "$meta.id"
    label 'process_very_high'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs:v0.4.1':
        'hkubal/clairs:v0.4.1' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(tumor_bam), path(tumor_bai),path(normal_bam), path(normal_bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_index)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.vcf.gz"), emit: somatic_vcf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def modelMap = [
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'hifi_revio'                        : 'hifi_revio_ss'
    
    ]
    
    // if method meta data is pb, default to the revio model else go to the map
    def model = modelMap.get(meta.basecall_model.toString().trim())

    if (!model in modelMap.keySet() ) {
        model = 'ont_r10_dorado_sup_5khz_ssrs'
        log.warn "Warning: ClairS-TO has no appropriate models for ${model} defaulting to dna_r10.4.1_e8.2_400bps_sup@v5.0.0 for Clair3"
    }
    else {
        log.info "Using ${model} model for Clair3"
    }
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
    /opt/bin/run_clairs \\
        --tumor_bam_fn $tumor_bam \\
        --normal_bam_fn $normal_bam \\
        --ref_fn $ref \\
        --threads ${task.cpus} \\
        --platform ${model} \\
        --output_dir . \\
        --use_heterozygous_indel_for_intermediate_phasing True \\
        --use_longphase_for_intermediate_haplotagging True \\
        --conda_prefix /opt/conda/envs/clairs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         clairs: \$( /opt/bin/run_clairs --version |& sed 's/ClairS v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         clairs: \$( /opt/bin/run_clairs --version |& sed 's/ClairS v//' )
    END_VERSIONS
    """
}
