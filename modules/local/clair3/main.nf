process CLAIR3 {
    tag "$meta.id"
    label 'process_very_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.2.0--py310h779eee5_0':
        'quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0' }"

    input:
    tuple val(meta) , path(bam), path(bai), path(model), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("*merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("*phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("*phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    tuple val("${task.process}"), val('clair3'), eval("run_clair3.sh  --version |& sed '1!d ; s/Clair3 v//'"), topic: versions, emit: versions_clair3

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    """
    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=$task.cpus \\
        --output=${prefix} \\
        --platform=$platform \\
        --model=$model \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz
    touch ${prefix}.phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merge_output.vcf.gz
    touch ${prefix}.merge_output.vcf.gz.tbi
    """
}
