process CLAIR3 {
    tag "$meta.id"
    label 'process_high'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clair3:v1.0.10':
        'hkubal/clair3:v1.0.10' }"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_index)


    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*merge_output.vcf.gz"), emit: germline_vcf

    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def modelMap = [
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'r1041_e82_400bps_sup_v500',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'r1041_e82_400bps_sup_v430',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'r1041_e82_400bps_sup_v420',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'r1041_e82_400bps_sup_v410',
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'r1041_e82_260bps_sup_v400',
        'hifi_revio'                        : 'hifi_revio'
    ]
    def model = modelMap.get(meta.basecall_model.toString().trim())
    def platform = (meta.platform == "pb") ? "hifi" : "ont"

    if (!model in modelMap.keySet() ) {
        model = 'r1041_e82_400bps_sup_v500'
        log.warn "Warning: ClairS-TO has no appropriate models for ${model} defaulting to dna_r10.4.1_e8.2_400bps_sup@v5.0.0 for Clair3"
    }
    else {
        log.info "Using ${model} model for Clair3"
    }
    
    // Download model command
    def download_prefix = ( model == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
    
    
    // Specify runtype for conda vs docker/singularity
    // Both the initial clair3 call and the model_path call will not work with conda
    //TODO:

                
    """
    wget ${download_prefix}/${model}.tar.gz
    tar -xvzf ${model}.tar.gz
    
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${task.cpus} \\
        --platform="${platform}" \\
        --model_path="${model}" \\
        --use_longphase_for_intermediate_phasing \\
        --output .               

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /opt/bin/run_clair3.sh --version |& sed 's/Clair3 v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merge_output.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /opt/bin/run_clair3.sh --version |& sed 's/Clair3 v//' )
    END_VERSIONS
    """
}
