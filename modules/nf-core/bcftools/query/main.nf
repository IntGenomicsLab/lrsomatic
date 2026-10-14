process BCFTOOLS_QUERY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'quay.io/biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path regions
    path targets
    path samples

    output:
    tuple val(meta), path("*.${suffix}.gz"), emit: output
    tuple val(meta), path("*.${suffix}.gz.tbi"), emit: index
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "txt"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file = samples ? "--samples-file ${samples}" : ""
    """
    bcftools query \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.${suffix}
    bgzip -c ${prefix}.${suffix} > ${prefix}.${suffix}.gz
    tabix -s 1 -b 2 -e 2 ${prefix}.${suffix}.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "txt"
    """
    touch ${prefix}.${suffix}
    bgzip -c ${prefix}.${suffix} > ${prefix}.${suffix}.gz
    touch ${prefix}.${suffix}.gz.tbi
    tabix -s 1 -b 2 -e 2 ${prefix}.${suffix}.gz
    """
}
