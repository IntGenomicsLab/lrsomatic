
process CLAIRSTO {
    tag "$meta.id"
    label 'process_very_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/hkubal/clairs-to:v0.4.2':
        'docker.io/hkubal/clairs-to:v0.4.2' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model), path(pon_vcfs), val(pon_flags)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("indel.vcf.gz"),      emit: indel_vcf
    tuple val(meta), path("indel.vcf.gz.tbi"),  emit: indel_tbi
    tuple val(meta), path("snv.vcf.gz"),        emit: snv_vcf
    tuple val(meta), path("snv.vcf.gz.tbi"),    emit: snv_tbi
    tuple val("${task.process}"), val('clairsto'), eval("run_clairs_to  --version |& sed '1!d ; s/run_clairs_to //'"), topic: versions, emit: versions_clairsto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def conda_prefix = workflow.containerEngine == 'singularity' ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''
    def pon_string   = pon_vcfs.join(',')
    def flags_string = pon_flags.join(',')

    """
    /opt/bin/run_clairs_to \
        --tumor_bam_fn $tumor_bam \\
        --ref_fn $reference \\
        --platform $model \\
        --threads $task.cpus \\
        --output_dir . \\
        --sample_name ${prefix} \\
        --panel_of_normals ${pon_string} \\
        --panel_of_normals_require_allele_matching ${flags_string} \\
        $conda_prefix \\
        $args \\
    """

    stub:
    """
    mkdir -p output
    echo "" | gzip > snv.vcf.gz
    touch snv.vcf.gz.tbi
    echo "" | gzip > indel.vcf.gz
    touch indel.vcf.gz.tbi
    """
}
