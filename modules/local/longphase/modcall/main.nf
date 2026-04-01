process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83fce1d397cf71705cc096fc0e0e52f7013bdd471ef68ee53ae765688e5c439c/data':
        'community.wave.seqera.io/library/longphase_samtools:8c61296cae7a5fc0' }"

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
