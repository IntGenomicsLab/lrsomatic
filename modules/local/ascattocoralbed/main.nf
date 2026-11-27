process ASCAT_TO_CORAL_BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12':
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(cnvs)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*_coral_cn.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fai_arg = fai ? "--fai ${fai}" : ''
    """
    ascat_to_coral_bed.py \\
        --cnvs ${cnvs} \\
        ${fai_arg} \\
        --output ${prefix}_coral_cn.bed \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_coral_cn.bed
    """
}
