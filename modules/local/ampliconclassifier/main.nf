process AMPLICONCLASSIFIER {
    tag "$meta.id"
    label 'process_medium'

    // No conda: bioconda's `ampliconclassifier` recipe is stuck at 0.4.14 (2023) and
    // predates CoRAL support; this image carries the CHM13 fork of v2.0.0. See meta.yml
    container "docker.io/robertaforsyth/ampliconclassifier:2.0.0-chm13-cdeaa63"

    input:
    tuple val(meta), path(reconstruction)
    tuple val(meta2), path(data_repo, stageAs: 'aa_data_repo')
    val(ac_ref)

    output:
    tuple val(meta), path("*_amplicon_classification_profiles.tsv"), emit: classification, optional: true
    tuple val(meta), path("*_gene_list.tsv"), emit: gene_list, optional: true
    tuple val(meta), path("*_ecDNA_counts.tsv"), emit: ecdna_counts, optional: true
    tuple val(meta), path("*_result_table.tsv"), emit: result_table, optional: true
    tuple val(meta), path("*_classification_bed_files", type: 'dir'), emit: bed_files, optional: true
    tuple val(meta), path("*_SV_summaries", type: 'dir'), emit: sv_summaries, optional: true
    tuple val(meta), path("*_annotated_cycles_files", type: 'dir'), emit: annotated_cycles, optional: true
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('ampliconclassifier'), eval("amplicon_classifier.py --version"), topic: versions, emit: versions_ampliconclassifier

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "AMPLICONCLASSIFIER does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export AA_DATA_REPO=\$(readlink -f aa_data_repo)

    amplicon_classifier.py \\
        --ref ${ac_ref} \\
        --AA_results ${reconstruction} \\
        -o ${prefix} \\
        ${args} \\
        > ${prefix}_classifier.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_amplicon_classification_profiles.tsv
    touch ${prefix}_gene_list.tsv
    touch ${prefix}_ecDNA_counts.tsv
    touch ${prefix}_result_table.tsv
    touch ${prefix}_classifier.log
    mkdir -p ${prefix}_classification_bed_files ${prefix}_SV_summaries ${prefix}_annotated_cycles_files
    """
}
