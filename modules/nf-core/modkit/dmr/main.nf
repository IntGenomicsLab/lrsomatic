process MODKIT_DMR {
    tag "${meta.id}"
    label 'process_medium'

    // Conda is not supported: matches MODKIT_PILEUP's own patch (the guard in `script:` stops
    // conda/mamba runs) -- kept on the same patched build as pileup for a single consistent
    // modkit binary across the whole DMR chain, rather than mixing stock 0.6.1 (dmr pair) with
    // the patched 0.6.4 (pileup). Revert to the biocontainer once the patch is upstreamed/released.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/modkit-sif:0.6.4-pacbiofix-6e0afa2'
        : 'ghcr.io/ljwharbers/modkit:0.6.4-pacbiofix-6e0afa2'}"

    input:
    tuple val(meta), path(bedmethyl_a), path(bedmethyl_a_tbi)
    tuple val(meta2), path(bedmethyl_b), path(bedmethyl_b_tbi)
    tuple val(meta3), path(regions_bed)
    tuple val(meta4), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MODKIT_DMR does not support Conda: the patched container is Docker/Singularity/Apptainer-only. Use one of those, or --skip_dmr."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-r ${regions_bed}" : ''
    // fai is never referenced directly -- modkit dmr pair has no --fai flag and looks for
    // <fasta>.fai next to --ref itself. Declaring it as an input only stages it alongside
    // fasta so modkit finds and reuses it instead of building one from scratch, which for a
    // multi-GB genome is a real, repeated cost across every task invocation. Confirmed by
    // direct testing: with no .fai present, modkit dmr pair creates one; with one already
    // staged next to fasta, it's read (unmodified mtime) rather than rebuilt.
    """
    modkit \\
        dmr pair \\
        -a ${bedmethyl_a} \\
        -b ${bedmethyl_b} \\
        ${regions} \\
        --ref ${fasta} \\
        -o ${prefix}.bed \\
        --log-filepath ${prefix}.log \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.log
    """
}
