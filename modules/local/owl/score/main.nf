process OWL_SCORE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/owl:0.4.0--h9ee0642_1' :
        'biocontainers/owl:0.4.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(profile)

    output:
    tuple val(meta), path("*.owl-scores.txt")        , emit: scores
    tuple val(meta), path("*.owl-motif-counts.txt")  , emit: motif_counts
    tuple val(meta), path("*.owl-phased-sites.txt")  , emit: phased_sites
    tuple val(meta), path("*.owl-phased-summary.txt"), emit: phased_summary
    tuple val("${task.process}"), val('owl'), eval("owl --version | sed 's/^owl //'"), topic: versions, emit: versions_owl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    owl score \\
        --file ${profile} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.owl-scores.txt
    touch ${prefix}.owl-motif-counts.txt
    touch ${prefix}.owl-phased-sites.txt
    touch ${prefix}.owl-phased-summary.txt
    """
}
