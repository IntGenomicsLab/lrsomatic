process CORAL_SEED {
    tag "$meta.id"
    label 'process_low'

    // No conda: CoRAL is not packaged on bioconda (the `coral` recipe there is an
    // unrelated RNA-seq tool), and this image carries the CHM13 fork. See meta.yml
    container "docker.io/robertaforsyth/coral:3.0.0-chm13-847f3d4"

    input:
    tuple val(meta), path(cn_seg), path(bam), path(bai)
    val(coral_ref)

    output:
    tuple val(meta), path("*_CNV_SEEDS.bed"), emit: seeds
    tuple val("${task.process}"), val('coral'), eval("python -c 'import importlib.metadata as m; print(m.version(\"CoRAL\"))'"), topic: versions, emit: versions_coral

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CORAL_SEED does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coral seed \\
        --cn-seg ${cn_seg} \\
        --ref ${coral_ref} \\
        --lr-bam ${bam} \\
        --output-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A non-empty seed, so the stub exercises the reconstruct path rather than the
    // empty-seed branch. Use --coral_gain to force the empty case in a test.
    """
    printf 'chr1\\t100000\\t400000\\t8\\n' > ${prefix}_CNV_SEEDS.bed
    """
}
