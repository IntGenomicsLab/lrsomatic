process CRAMINO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cramino:1.3.0--h3dc2dae_0':
        'biocontainers/cramino:1.3.0--h3dc2dae_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('cramino'), eval("cramino --version |& sed '1!d ; s/cramino //'"), topic: versions, emit: versions_cramino

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cramino $args $bam > ${prefix}_cramino.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cramino.txt
    """
}
