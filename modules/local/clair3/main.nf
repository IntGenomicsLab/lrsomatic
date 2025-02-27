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

process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clair3:v1.0.10':
        'hkubal/clair3:v1.0.10' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam), path(bam_bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_index)


    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("merge_output.vcf.gz"), emit: germline_vcf

    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def modelMap = [
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'r1041_e82_400bps_sup_v500',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'r1041_e82_400bps_sup_v430',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'r1041_e82_400bps_sup_v420',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'r1041_e82_400bps_sup_v410',
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'r1041_e82_260bps_sup_v400',
        'hifi_revio'                        : 'hifi_revio'
    ]
    def model = modelMap.get(meta.basecall_model.toString().trim())
    def platform = (meta.platform == "pb") ? "hifi" : "ont"

    if (!model in modelMap.keySet() ) {
        model = 'r1041_e82_400bps_sup_v500'
        log.warn "Warning: ClairS-TO has no appropriate models for ${model} defaulting to dna_r10.4.1_e8.2_400bps_sup@v5.0.0 for Clair3"
    }
    else {
        log.info "Using ${model} model for Clair3"
    }
    
    // Download model command
    def download_prefix = ( model == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
    
    
    // Specify runtype for conda vs docker/singularity
    // Both the initial clair3 call and the model_path call will not work with conda
    //TODO:

                
    """
    wget ${download_prefix}/${model}.tar.gz
    tar -xvzf ${model}.tar.gz
    
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${task.cpus} \\
        --platform="${platform}" \\
        --model_path="${model}" \\
        --use_longphase_for_intermediate_phasing \\
        --output .               

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /opt/bin/run_clair3.sh --version |& sed 's/Clair3 v//' )
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
        clair3: \$( /opt/bin/run_clair3.sh --version |& sed 's/Clair3 v//' )
    END_VERSIONS
    """
}
