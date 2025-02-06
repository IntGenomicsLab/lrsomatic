process CRAMINO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cramino:0.16.0--h3dc2dae_0':
        'biocontainers/cramino:0.16.0--h3dc2dae_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    




    """
    cramino $args $bam > ${prefix}_cramino.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version |& sed '1!d ; s/cramino //')
    END_VERSIONS
    """





    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cramino.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version |& sed '1!d ; s/cramino //')
    END_VERSIONS
    """
}
