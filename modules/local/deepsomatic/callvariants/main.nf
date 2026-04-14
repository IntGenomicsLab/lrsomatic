process DEEPSOMATIC_CALLVARIANTS {
    tag "$meta.id"
    label "${params.use_gpu ? 'process_gpu_high' : 'process_very_high'}"

    //Conda is not supported at the moment
    container params.use_gpu ? "docker.io/google/deepsomatic:1.7.0-gpu" : "docker.io/google/deepsomatic:1.7.0"

    input:
    tuple val(meta), path(make_examples_tfrecords)

    output:
    tuple val(meta), path("${prefix}.call-*-of-*.tfrecord.gz")    , emit: call_variants_tfrecords
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

    def matcher = make_examples_tfrecords[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    if (!matcher.matches()) {
        throw new IllegalArgumentException("tfrecord baseName '" + make_examples_tfrecords[0].baseName + "' doesn't match the expected pattern")
    }
    def examples_tfrecord_name = matcher[0][1]
    def shardCount = matcher[0][2]
    // Reconstruct the logical name - ${tfrecord_name}@<shard_count>.gz
    def examples_tfrecords_logical_name = "${examples_tfrecord_name}@${shardCount}.gz"

    """
    /opt/deepvariant/bin/call_variants \\
        ${args} \\
        --outfile "${prefix}.call.tfrecord.gz" \\
        --examples "${examples_tfrecords_logical_name}"

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.call-00000-of-00001.tfrecord.gz

    """
}
