process CORAL_CYCLE {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'ignore' }

    container "docker.io/robertaforsyth/coral:3.0.0-chm13-847f3d4"

    input:
    tuple val(meta), path(reconstruction)

    output:
    tuple val(meta), path("cycles"), emit: reconstruction
    tuple val(meta), path("cycles/*_amplicon*_cycles.txt"), emit: cycles, optional: true
    tuple val("${task.process}"), val('coral'), eval("python -c 'import importlib.metadata as m; print(m.version(\"CoRAL\"))'"), topic: versions, emit: versions_coral

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CORAL_CYCLE does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Graphs are copied alongside the re-extracted cycles so AmpliconClassifier
    // still finds a matching graph/cycles/summary set in one directory.
    """
    mkdir -p cycles
    cp ${reconstruction}/*_graph.txt ${reconstruction}/*_summary.txt cycles/ 2>/dev/null || true

    coral cycle_all \\
        --bp-dir ${reconstruction} \\
        --output-prefix cycles/${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p cycles
    touch cycles/${prefix}_amplicon1_cycles.txt
    touch cycles/${prefix}_amplicon1_graph.txt
    touch cycles/${prefix}_summary.txt
    """
}
