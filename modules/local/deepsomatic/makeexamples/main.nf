process DEEPSOMATIC_MAKEEXAMPLES {
    tag "$meta.id"
    label 'process_very_high'
    label 'process_long'

    //Conda is not supported at the moment
    container params.use_gpu ? "docker.io/google/deepsomatic:1.7.0-gpu" : "docker.io/google/deepsomatic:1.7.0"

    input:
    tuple val(meta), path(normal_input), path(normal_index), path(tumor_input), path(tumor_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(ds_pon)

    output:
    tuple val(meta), path("${prefix}.examples.tfrecord-*-of-*.gz{,.example_info.json}")         , emit: examples
    tuple val(meta), path("${prefix}.gvcf.tfrecord-*-of-*.gz")                                  , emit: gvcf, optional:true
    tuple val(meta), path("${prefix}_call_variant_outputs.tfrecord-*-of-*.gz", arity: "0..*")   , emit: small_model_calls
    tuple val("${task.process}"), val('deepsomatic'), val('1.7.0'), topic: versions, emit: versions_deepsomatic

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def normalReadsArg = (normal_input?.toString() && normal_input.toString() != '[]') ? "--reads_normal \"${normal_input}\"" : ""
    def normalSampleArg = (normal_input?.toString() && normal_input.toString() != '[]') ? "--sample_name_normal \"${prefix}_normal\"" : ""
    def gvcf_arg = params.generate_gvcf ? "--gvcf \"./${prefix}.gvcf.tfrecord@${task.cpus}.gz\"" : ""
    def isTumorOnly = !(normal_input?.toString() && normal_input.toString() != '[]')
    def ponArg = ""
    if (ds_pon?.toString() && ds_pon.toString() != '[]') {
        // User-supplied PON: staged into work dir; works in both tumor-only and paired modes
        def ponPaths = ds_pon instanceof List
            ? ds_pon.findAll { !it.toString().endsWith('.tbi') }.collect { "\"${it}\"" }.join(',')
            : "\"${ds_pon}\""
        ponArg = "--population_vcfs ${ponPaths}"
    } else if (isTumorOnly) {
        // No user PON in tumor-only mode: fall back to container-bundled defaults
        ponArg = '--population_vcfs "/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz","AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz","PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz","PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz"'
    }
    // In paired mode with no user PON, ponArg stays "" (no --population_vcfs, matching prior behaviour)

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples_somatic \\
        --mode calling \\
        --ref "${fasta}" \\
        --reads_tumor "${tumor_input}" \\
        ${normalReadsArg} \\
        --sample_name_tumor "${prefix}" \\
        ${normalSampleArg} \\
        --examples "./${prefix}.examples.tfrecord@${task.cpus}.gz" \\
        ${gvcf_arg} \\
        ${ponArg} \\
        ${args} \\
        --task {}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf -v SHARD_COUNT "%04d" ${task.cpus}
    for i in \$( seq -f "%04g" 0 ${task.cpus-1} )
    do
        echo "" | gzip > ${prefix}.examples.tfrecord-\$i-of-\$SHARD_COUNT.gz
        touch ${prefix}.examples.tfrecord-\$i-of-\$SHARD_COUNT.gz.example_info.json
        echo "" | gzip > ${prefix}.gvcf.tfrecord-\$i-of-\$SHARD_COUNT.gz
    done
    """
}
