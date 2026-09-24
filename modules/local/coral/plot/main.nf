process CORAL_PLOT {
    tag "$meta.id"
    label 'process_medium'

    // Plotting is cosmetic: never let it fail a run that reconstructed successfully
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'ignore' }

    container "docker.io/robertaforsyth/coral:3.0.0-chm13-847f3d4"

    input:
    tuple val(meta), path(reconstruction), path(bam), path(bai)
    val(coral_ref)

    output:
    tuple val(meta), path("*_amplicon*.{pdf,png}"), emit: plots, optional: true
    tuple val("${task.process}"), val('coral'), eval("python -c 'import importlib.metadata as m; print(m.version(\"CoRAL\"))'"), topic: versions, emit: versions_coral

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CORAL_PLOT does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coral plot_all \\
        --ref ${coral_ref} \\
        --reconstruction-dir ${reconstruction} \\
        --bam ${bam} \\
        --output-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_amplicon1_graph.png
    touch ${prefix}_amplicon1_cycles.png
    """
}
