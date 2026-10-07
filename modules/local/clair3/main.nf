process CLAIR3 {
    tag "$meta.id"
    label "${params.use_gpu ? 'process_gpu_very_high' : 'process_very_high'}"

    conda "${moduleDir}/environment.yml"
    // No GPU image was published for 2.0.3; 2.0.2 differs only by the --gender option, which is not used here
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (params.use_gpu ? 'docker://docker.io/hkubal/clair3:v2.0.2_gpu' : 'https://depot.galaxyproject.org/singularity/clair3:2.0.3--py311hbc58adc_0') :
        (params.use_gpu ? 'docker.io/hkubal/clair3:v2.0.2_gpu' : 'quay.io/biocontainers/clair3:2.0.3--py311hbc58adc_0') }"

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
    // An empty model input means the model is bundled with the image: models/ next to run_clair3.sh
    // in the biocontainer, /opt/models in the HKU Docker image
    def model_path = model ? "${model}" : "\${CLAIR3_MODELS}/${meta.clair3_model}"

    """
    ${use_gpu ? 'export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}' : ':'}
    CLAIR3_MODELS=\$(dirname \$(command -v run_clair3.sh))/models
    [ -d "\${CLAIR3_MODELS}" ] || CLAIR3_MODELS=/opt/models

    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model_path=${model_path} \\
        --sample_name=${prefix} \\
        ${use_gpu ? '--use_gpu --device=cuda:0' : ''} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz
    touch ${prefix}.phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merge_output.vcf.gz
    touch ${prefix}.merge_output.vcf.gz.tbi
    """
}
