process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_very_high'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf")                       , optional: true, emit: paf
    tuple val(meta), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta), path("*.cram")                      , optional: true, emit: cram
    tuple val(meta), path("*.{bam,cram}.${bam_index_extension}"), optional: true, emit: index
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // bam_format 'cram' writes a reference-compressed CRAM instead of a BAM
    def out_ext = bam_format == 'cram' ? 'cram' : 'bam'
    def cram_opts = bam_format == 'cram' ? "-O cram --reference ${reference}" : ''
    def bam_index = bam_index_extension ? "${prefix}.${out_ext}##idx##${prefix}.${out_ext}.${bam_index_extension} --write-index" : "${prefix}.${out_ext}"
    def bam_output = bam_format ? "-a | samtools sort -@ ${task.cpus-1} ${cram_opts} -o ${bam_index} ${args2}" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    // A CRAM input is decoded against `reference` only if no contig it shares with it has a different M5. Given a
    // mismatched --reference, an embedded-reference CRAM decodes to wrong bases with exit 0; without it, htslib uses the
    // embedded reference, the @SQ UR path or REF_PATH. Unaligned CRAM has no @SQ lines and needs no reference.
    def cram_input = "${reads.extension}" == 'cram' && reference
    def cram_ref_check = cram_input ? """reset_ref="--reference ${reference}"
    samtools view -H ${reads} | grep '^@SQ' > input_sq.txt || true
    if [ -s input_sq.txt ] && samtools dict ${reference} | grep '^@SQ' | awk -F'\\t' '
        function tag(t,   i) { for (i = 2; i <= NF; i++) if (substr(\$i, 1, 3) == t ":") return substr(\$i, 4); return "" }
        FNR == NR { m5[tag("SN")] = tag("M5"); next }
        (tag("SN") in m5) && tag("M5") != "" && tag("M5") != m5[tag("SN")] { bad = 1 }
        END { exit !bad }' - input_sq.txt; then
        echo "WARNING: ${reads} is aligned to a different assembly than ${reference}; decoding it with its embedded reference, @SQ UR or REF_PATH instead" >&2
        reset_ref=""
    fi
    """ : ''
    // samtools reset (1.23) takes no --reference, so a CRAM is decoded by samtools view first
    def samtools_reset = cram_input
        ? "samtools view -u --threads ${task.cpus-1} \$reset_ref $reads | samtools reset --threads ${task.cpus-1} $args3 -"
        : "samtools reset --threads ${task.cpus-1} $args3 $reads"
    def samtools_reset_fastq = bam_input ? "$samtools_reset | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("Error: minimap2/align BAM input mode requires reference") : reads)
    """
    ${cram_ref_check}$samtools_reset_fastq \\
    minimap2 \\
        ${args} \\
        -t ${task.cpus} \\
        ${target} \\
        ${query} \\
        ${cigar_paf} \\
        ${set_cigar_bam} \\
        ${bam_output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_ext = bam_format == 'cram' ? 'cram' : 'bam'
    def output_file = bam_format ? "${prefix}.${out_ext}" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.${out_ext}.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    if(bam_input && !reference) {
        error("Error: minimap2/align BAM input mode requires reference!")
	}
    """
    touch $output_file
    ${bam_index}
    """
}
