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
    def isTumorOnly = !(meta.paired_data)

    // Build list of PON VCF file paths (excluding .tbi index files)
    def ponFiles = []
    if (ds_pon?.toString() && ds_pon.toString() != '[]') {
        ponFiles = (ds_pon instanceof List)
            ? ds_pon.findAll { !it.toString().endsWith('.tbi') }
            : [ds_pon]
    }
    def nPonFiles = ponFiles.size()
    def ponArrayLiteral = ponFiles.collect { "${it}" }.join(' ')

    // Shell block to prepare the PON VCF (merge if multiple, copy if single, skip if none)
    // Runs before make_examples_somatic so the result is available as merged_pon.vcf.gz
    def ponPrepareBlock = (isTumorOnly && nPonFiles > 0) ? """
# Prepare PON VCF for --population_vcfs: merge multiple databases into one sorted+indexed file,
# or copy a single VCF. DeepSomatic requires no chromosome overlap across population VCFs.
_PON_VCFS=( ${ponArrayLiteral} )
if [ \${#_PON_VCFS[@]} -gt 1 ]; then
    gzip -dc "\${_PON_VCFS[0]}" | grep '^##fileformat' > _pon_hdr.txt
    for vcf in "\${_PON_VCFS[@]}"; do gzip -dc "\$vcf" | grep '^##' | grep -v '^##fileformat'; done | sort -u >> _pon_hdr.txt
    gzip -dc "\${_PON_VCFS[0]}" | grep '^#CHROM' >> _pon_hdr.txt
    for vcf in "\${_PON_VCFS[@]}"; do gzip -dc "\$vcf" | grep -v '^#'; done \\
        | sort -t\$'\\t' -k1,1V -k2,2n | uniq > _pon_data.txt
    cat _pon_hdr.txt _pon_data.txt | bgzip -c > merged_pon.vcf.gz
    rm _pon_hdr.txt _pon_data.txt
else
    cp "\${_PON_VCFS[0]}" merged_pon.vcf.gz
fi
tabix -p vcf merged_pon.vcf.gz
""" : ""

    // --population_vcfs argument for make_examples_somatic
    def ponArg = ""
    if (isTumorOnly) {
        ponArg = nPonFiles > 0
            ? '--population_vcfs "merged_pon.vcf.gz"'
            : '--population_vcfs "/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz "'
    }
    // In paired mode ponArg stays "" (no --population_vcfs, matching prior behaviour)

    """
    ${ponPrepareBlock}
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
