process DEEPSOMATIC_POSTPROCESSVARIANTS {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepsomatic:1.7.0"

    input:
    tuple val(meta), path(variant_calls_tfrecord_files), path(gvcf_tfrecords), val(small_model_calls), val(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(pon_vcf)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),                                        emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.{tbi,csi}"),                              emit: vcf_index
    tuple val(meta), path("${prefix}.g.vcf.gz"),                                      emit: gvcf,                  optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.{tbi,csi}"),                            emit: gvcf_index,            optional: true
    tuple val("${task.process}"), val('deepsomatic'), val('1.7.0'), topic: versions,  emit: versions_deepsomatic

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def regions = intervals ? "--regions ${intervals}" : ""
    def variant_calls_tfrecord_name = variant_calls_tfrecord_files[0].name.replaceFirst(/-\d{5}-of-\d{5}/, "")

    def isTumorOnly = !(meta.paired_data)


    if (gvcf_tfrecords) {
        def gvcf_matcher = gvcf_tfrecords[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
        if (!gvcf_matcher.matches()) {
            throw new IllegalArgumentException("tfrecord baseName '" + gvcf_tfrecords[0].baseName + "' doesn't match the expected pattern")
        }
    }

    // The following block determines whether the small model was used, and if so, adds the variant calls from it
    // to the argument --small_model_cvo_records.
    def small_model_arg = ""
    if (small_model_calls && small_model_calls.size() > 0) {
        def small_model_matcher = (small_model_calls[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/)
        if (!small_model_matcher.matches()) {
            throw new IllegalArgumentException("tfrecord baseName '" + small_model_calls[0].baseName + "' doesn't match the expected pattern")
        }
        def small_model_tfrecord_name = small_model_matcher[0][1]
        def small_model_shardCount = small_model_matcher[0][2]
        // Reconstruct the logical name. Example: test_call_variant_outputs.examples.tfrecord@12.gz
        def small_model_tfrecords_logical_name = "${small_model_tfrecord_name}@${small_model_shardCount}.gz"
        small_model_arg = "--small_model_cvo_records ${small_model_tfrecords_logical_name}"
    }

    // Build list of PON VCF file paths (excluding .tbi index files)
    def ponFiles = []
    if (pon_vcf?.toString() && pon_vcf.toString() != '[]') {
        ponFiles = (pon_vcf instanceof List)
            ? pon_vcf.findAll { f -> !f.toString().endsWith('.tbi') }
            : [pon_vcf]
    }
    def nPonFiles = ponFiles.size()
    def ponArrayLiteral = ponFiles.collect { f -> "${f}" }.join(' ')

    // Shell block to prepare the PON VCF for --pon_filtering (merge if multiple, copy if single)
    def ponPrepareBlock = (isTumorOnly && nPonFiles > 0) ? """
# Prepare PON VCF for --pon_filtering: merge multiple databases into one sorted+indexed file,
# or copy a single VCF. DeepSomatic requires a single VCF for --pon_filtering.
_PON_VCFS=( ${ponArrayLiteral} )
if [ \${#_PON_VCFS[@]} -gt 1 ]; then
    gzip -dc "\${_PON_VCFS[0]}" | grep '^##fileformat' > _pon_hdr.txt
    for vcf in "\${_PON_VCFS[@]}"; do gzip -dc "\$vcf" | grep '^##' | grep -v '^##fileformat'; done | sort -T . -u >> _pon_hdr.txt
    gzip -dc "\${_PON_VCFS[0]}" | grep '^#CHROM' >> _pon_hdr.txt
    for vcf in "\${_PON_VCFS[@]}"; do gzip -dc "\$vcf" | grep -v '^#'; done \\
        | sort -T . -t\$'\\t' -k1,1V -k2,2n | uniq > _pon_data.txt
    cat _pon_hdr.txt _pon_data.txt | bgzip -c > merged_pon.vcf.gz
    rm _pon_hdr.txt _pon_data.txt
else
    cp "\${_PON_VCFS[0]}" merged_pon.vcf.gz
fi
tabix -p vcf merged_pon.vcf.gz
""" : ""

    // --pon_filtering argument for postprocess_variants (tumor-only only)
    def ponFilterArg = ""
    if (isTumorOnly) {
        ponFilterArg = nPonFiles > 0
            ? '--pon_filtering "merged_pon.vcf.gz"'
            : '--pon_filtering "/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz"'
    }
    // Paired samples: ponFilterArg stays "" (no PON filtering)

    """
    ${ponPrepareBlock}
    /opt/deepvariant/bin/postprocess_variants \\
        ${args} \\
        --ref "${fasta}" \\
        --infile "${variant_calls_tfrecord_name}" \\
        --outfile "${prefix}.vcf.gz" \\
        --process_somatic=true \\
        ${regions} \\
        ${small_model_arg} \\
        ${ponFilterArg} \\
        --cpus ${task.cpus}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    """
}
