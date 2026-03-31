process CLAIR3 {
    tag "$meta.id"
    label "${params.use_gpu ? 'process_gpu_very_high' : 'process_very_high'}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.use_gpu ? 'docker://hkubal/clair3-gpu:v1.2.0' : 'https://depot.galaxyproject.org/singularity/clair3:1.2.0--py310h779eee5_0') :
        (params.use_gpu ? 'docker.io/hkubal/clair3-gpu:v1.2.0' : 'quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0') }"

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
    prefix = task.ext.prefix ?: "${meta.id}"
    def use_gpu = task.ext.use_gpu as boolean

    """
    ${use_gpu ? 'export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}' : ':'}
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model=${model} \\
        --sample_name=${prefix} \\
        ${use_gpu ? '--use_gpu --device=cuda:0' : ''} \\
        ${args}
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
