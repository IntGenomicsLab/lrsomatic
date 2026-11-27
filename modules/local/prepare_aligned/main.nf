process PREPARE_ALIGNED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta), path(fai)
    val out_format   // 'bam' or 'cram'

    output:
    tuple val(meta), path("${prefix}.${out_ext}")           , emit: alignment
    tuple val(meta), path("${prefix}.${out_ext}.${out_idx}"), emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''    // samtools view: tags to drop (-x)
    def args2      = task.ext.args2 ?: ''   // samtools addreplacerg: output options
    def read_group = task.ext.read_group ?: error("PREPARE_ALIGNED needs ext.read_group (an @RG line with \\t separators)")
    prefix  = task.ext.prefix ?: "${meta.id}"
    out_ext = out_format == 'cram' ? 'cram' : 'bam'
    out_idx = out_format == 'cram' ? 'crai' : 'bai'
    if ("${reads}" == "${prefix}.${out_ext}") {
        error("PREPARE_ALIGNED: input and output are both named ${reads}; set a different ext.prefix")
    }
    // A CRAM is decoded with --fasta only if no contig it shares with it has a different M5 (see MINIMAP2_ALIGN)
    def cram_ref_check = reads.extension == 'cram' ? """
    if samtools dict ${fasta} | grep '^@SQ' | awk -F'\\t' '
        function tag(t,   i) { for (i = 2; i <= NF; i++) if (substr(\$i, 1, 3) == t ":") return substr(\$i, 4); return "" }
        FNR == NR { m5[tag("SN")] = tag("M5"); next }
        (tag("SN") in m5) && tag("M5") != "" && tag("M5") != m5[tag("SN")] { bad = 1 }
        END { exit bad }' - input_sq.txt; then
        decode_ref="--reference ${fasta}"
    fi
    """ : ''
    """
    samtools view -H ${reads} > input_header.sam
    grep '^@SQ' input_header.sam > input_sq.txt || true

    if [ ! -s input_sq.txt ]; then
        echo "ERROR: ${reads} has no @SQ lines, so it is not aligned. --skip_alignment needs reads aligned to --fasta; run without it to align them." >&2
        exit 1
    fi

    # Every contig of the input must be in --fasta with the same length
    if ! awk -F'\\t' '
        FNR == NR { len[\$1] = \$2; next }
        {
            sn = ""; ln = ""
            for (i = 2; i <= NF; i++) {
                if (substr(\$i, 1, 3) == "SN:") sn = substr(\$i, 4)
                if (substr(\$i, 1, 3) == "LN:") ln = substr(\$i, 4)
            }
            if (!(sn in len)) { print "  " sn ": not in --fasta"; bad = 1 }
            else if (len[sn] != ln) { print "  " sn ": length " ln " in the input, " len[sn] " in --fasta"; bad = 1 }
        }
        END { exit bad }' ${fai} input_sq.txt > sq_mismatch.txt; then
        echo "ERROR: ${reads} was not aligned to --fasta (${fasta}). Contigs that differ:" >&2
        head -n 20 sq_mismatch.txt >&2
        echo "Run without --skip_alignment to realign it, or use the --fasta it was aligned to." >&2
        exit 1
    fi

    if ! grep -q '^@HD.*SO:coordinate' input_header.sam; then
        echo "ERROR: ${reads} is not coordinate-sorted (no SO:coordinate in @HD). Sort it with samtools sort first." >&2
        exit 1
    fi

    decode_ref=""
    ${cram_ref_check}
    # The pipeline's read group replaces every input one, as MINIMAP2_ALIGN's -R does
    rg_line=\$(printf '%b' '${read_group}')
    rg_id=\$(printf '%s\\n' "\$rg_line" | tr '\\t' '\\n' | sed -n 's/^ID://p')
    { grep -v '^@RG' input_header.sam; printf '%s\\n' "\$rg_line"; } > new_header.sam

    samtools view --no-PG -u --threads ${task.cpus - 1} \$decode_ref ${args} ${reads} \\
        | samtools reheader -P new_header.sam - \\
        | samtools addreplacerg \\
            --threads ${task.cpus - 1} \\
            -R "\$rg_id" \\
            -m overwrite_all \\
            -O ${out_ext} \\
            --reference ${fasta} \\
            ${args2} \\
            -o ${prefix}.${out_ext} \\
            -

    samtools index -@ ${task.cpus - 1} ${prefix}.${out_ext}
    """

    stub:
    prefix  = task.ext.prefix ?: "${meta.id}"
    out_ext = out_format == 'cram' ? 'cram' : 'bam'
    out_idx = out_format == 'cram' ? 'crai' : 'bai'
    """
    touch ${prefix}.${out_ext}
    touch ${prefix}.${out_ext}.${out_idx}
    """
}
