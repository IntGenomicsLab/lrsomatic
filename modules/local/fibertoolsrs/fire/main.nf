process FIBERTOOLSRS_FIRE {
    tag "$meta.id"
    label 'process_very_high'
    label "${params.use_gpu ? 'process_gpu_very_high_memory' : 'process_high_memory'}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fibertools-rs:0.8.1--h3b373d1_0':
        'biocontainers/fibertools-rs:0.8.1--h3b373d1_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fibertoolsrs'), eval("ft --version |& sed 's/^fibertools-rs v//; s/\tgit-details.*//'"), topic: versions, emit: versions_fibertoolsrs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ft \\
        fire \\
        $args \\
        -t $task.cpus \\
        $bam \\
        ${prefix}_fire.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fire.bam
    """
}
