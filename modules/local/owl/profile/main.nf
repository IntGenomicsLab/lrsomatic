process OWL_PROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/owl:0.4.0--h9ee0642_1' :
        'biocontainers/owl:0.4.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(markers)

    output:
    tuple val(meta), path("*.owl-profile.txt"), emit: profile
    tuple val("${task.process}"), val('owl'), eval("owl --version | sed 's/^owl //'"), topic: versions, emit: versions_owl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // minimap2 writes no pbmm2 mg tag, so owl warns once per read; drop those lines from the log
    """
    status=0
    owl profile \\
        --bam ${bam} \\
        --regions ${markers} \\
        --sample ${prefix} \\
        ${args} \\
        > ${prefix}.owl-profile.txt \\
        2> owl-profile.log || status=\$?
    grep -v 'keeping read without MG' owl-profile.log >&2 || true
    [ \$status -eq 0 ]
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.owl-profile.txt
    """
}
