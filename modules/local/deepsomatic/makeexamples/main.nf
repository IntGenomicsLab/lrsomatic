process DEEPSOMATIC_MAKEEXAMPLES {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "docker.io/google/deepsomatic:1.7.0"

    input:
    tuple val(meta), path(normal_input), path(normal_index), path(tumor_input), path(tumor_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.examples.tfrecord-*-of-*.gz{,.example_info.json}"), emit: examples
    tuple val(meta), path("${prefix}.gvcf.tfrecord-*-of-*.gz"), emit: gvcf
    tuple val(meta), path("${prefix}_call_variant_outputs.tfrecord-*-of-*.gz", arity: "0..*"), emit: small_model_calls
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

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples_somatic \\
        --mode calling \\
        --ref "${fasta}" \\
        --reads_tumor "${tumor_input}" \\
        ${normalReadsArg} \\
        --sample_name_tumor "${prefix}_tumor" \\
        ${normalSampleArg} \\
        --examples "./${prefix}.examples.tfrecord@${task.cpus}.gz" \\
        --gvcf "./${prefix}.gvcf.tfrecord@${task.cpus}.gz" \\
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
