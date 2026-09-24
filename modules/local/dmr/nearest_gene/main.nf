process DMR_NEAREST_GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(dmr_bed)
    path gencode_gene_bed

    output:
    tuple val(meta), path("*.nearestGene.bed"), emit: nearest_gene_bed
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-D a -t first'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sort -k1,1 -k2,2n ${dmr_bed} > dmr.sorted.bed

    bedtools closest -a dmr.sorted.bed -b ${gencode_gene_bed} ${args} > ${prefix}.nearestGene.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nearestGene.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
