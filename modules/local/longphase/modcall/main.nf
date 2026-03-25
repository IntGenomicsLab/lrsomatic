process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_very_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0184a9a36d8612fbae38bbaad7b52f03b815ad17673740e107cf1f267a1f15d/data':
        'community.wave.seqera.io/library/htslib_longphase:3071e61356fc25a4' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.vcf")       , emit: mod_vcf
    tuple val(meta), path("*.log")       , emit: log , optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    longphase \\
        modcall \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        -b ${bam} \\
        --out-prefix ${prefix}

    if [ -f "${prefix}.out" ]; then
        mv ${prefix}.out ${prefix}.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def log = args.contains('--log') ? "touch ${prefix}.log" : ''
    """
    touch ${prefix}.vcf
    ${log}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """
}
