process FIBERTOOLSRS_QC {
    tag "$meta.id"
    label 'process_very_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fibertools-rs:0.8.1--h3b373d1_0':
        'biocontainers/fibertools-rs:0.8.1--h3b373d1_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: qc_txt
    tuple val("${task.process}"), val('fibertoolsrs'), eval("ft --version |& sed '1!d ; s/ft //'"), topic: versions, emit: versions_fibertoolsrs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ft \\
        qc \\
        $args \\
        -t $task.cpus \\
        $bam \\
        ${prefix}_qc.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_qc.txt
    """
}
