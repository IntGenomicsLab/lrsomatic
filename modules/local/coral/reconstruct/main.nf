process CORAL_RECONSTRUCT {
    tag "$meta.id"
    label 'process_high'

    // A single unsolvable amplicon should not fail a whole cohort; the subworkflow
    // warns on the missing output rather than letting the report slot go silent.
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'ignore' }

    container "docker.io/robertaforsyth/coral:3.0.0-chm13-847f3d4"

    input:
    tuple val(meta), path(seeds), path(cn_seg), path(bam), path(bai)

    output:
    tuple val(meta), path("reconstruct"), emit: reconstruction
    tuple val(meta), path("reconstruct/*_amplicon*_graph.txt"), emit: graphs, optional: true
    tuple val(meta), path("reconstruct/*_amplicon*_cycles.txt"), emit: cycles, optional: true
    tuple val(meta), path("reconstruct/*_summary.txt"), emit: summary, optional: true
    tuple val(meta), path("reconstruct/*_reconstruct.log"), emit: log, optional: true
    tuple val("${task.process}"), val('coral'), eval("python -c 'import importlib.metadata as m; print(m.version(\"CoRAL\"))'"), topic: versions, emit: versions_coral

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CORAL_RECONSTRUCT does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // CoRAL derives its output directory by splitting --output-prefix on '/', and
    // AmpliconClassifier pairs graph/cycles/summary by prefix within one directory.
    """
    mkdir -p reconstruct

    coral reconstruct \\
        --lr-bam ${bam} \\
        --cnv-seed ${seeds} \\
        --cn-seg ${cn_seg} \\
        --output-prefix reconstruct/${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p reconstruct
    touch reconstruct/${prefix}_amplicon1_graph.txt
    touch reconstruct/${prefix}_amplicon1_cycles.txt
    touch reconstruct/${prefix}_summary.txt
    touch reconstruct/${prefix}_reconstruct.log
    """
}
