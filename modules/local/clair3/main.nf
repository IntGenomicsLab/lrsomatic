process CLAIR3 {
    tag "$meta.id"
    label 'process_very_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clair3:v1.0.10':
        'hkubal/clair3:v1.0.10' }"

    input:
    tuple val(meta), path(bam), path(bai), val(packaged_model), path(user_model), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("*merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("*phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("*phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = packaged_model
    //if (!user_model) {
    //    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //        model = "\${CONDA_PREFIX}/bin/models/${packaged_model}"
    //    }
    //    else {
    //        model = "/usr/local/bin/models/$packaged_model"
    //    }
    //}
    //if (!packaged_model) {
    //    model = "$user_model"
    //}
    //if (packaged_model && user_model) {
    //    error "Two models specified $user_model and $packaged_model, specify one of them."
    //}
    // TODO: fix the channel structure so you don't have to do this
    def download_prefix = ( model == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wget ${download_prefix}/${packaged_model}.tar.gz
    tar -xvzf ${model}.tar.gz

    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=$task.cpus \\
        --output=. \\
        --platform=$platform \\
        --model=$model \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh  --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz
    touch ${prefix}.phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merge_output.vcf.gz
    touch ${prefix}.merge_output.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """
}
