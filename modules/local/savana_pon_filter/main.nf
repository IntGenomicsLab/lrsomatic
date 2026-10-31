process SAVANA_PON_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12':
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(pon)

    output:
    tuple val(meta), path("${prefix}.classified.somatic.vcf"), emit: vcf
    tuple val(meta), path("${prefix}.pon_filter_summary.tsv") , emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // The output keeps the .classified.somatic.vcf suffix: the ReConPlot wrapper finds SAVANA's SV file by it
    prefix = task.ext.prefix ?: "${meta.id}.pon_filtered"
    """
    savana_pon_filter.py \\
        --vcf ${vcf} \\
        --pon ${pon} \\
        --output ${prefix}.classified.somatic.vcf \\
        --summary ${prefix}.pon_filter_summary.tsv \\
        --sample ${meta.id}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.pon_filtered"
    """
    touch ${prefix}.classified.somatic.vcf
    touch ${prefix}.pon_filter_summary.tsv
    """
}
