process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // 2.0.2, not 2.0.1: 2.0.1 read its region iterator with sam_itr_multi_next, which on CRAM input returns
    // fewer reads, differently every run (~28 % fewer sites genome-wide on a PacBio pair); 2.0.2 reads
    // CRAM like BAM and is deterministic. The script needs only longphase, so the plain biocontainer
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longphase:2.0.2--h4e109e1_0':
        'quay.io/biocontainers/longphase:2.0.2--h4e109e1_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.vcf")       , emit: mod_vcf
    tuple val(meta), path("*.log")       , emit: log , optional: true
    tuple val("${task.process}"), val('longphase'), eval("longphase --version | head -n 1 | sed 's/Version: //'"), topic: versions, emit: versions_longphase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    longphase \\
        modcall \\
        $args \\
        --threads 1 \\
        -o ${prefix} \\
        --reference ${fasta} \\
        -b ${bam} \\
        --out-prefix ${prefix}

    if [ -f "${prefix}.out" ]; then
        mv ${prefix}.out ${prefix}.log
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def log = args.contains('--log') ? "touch ${prefix}.log" : ''
    """
    touch ${prefix}.vcf
    ${log}
    """
}
