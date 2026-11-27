process METAEXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)   // the pipeline reference for CRAM input, [[:], []] otherwise

    output:
    tuple val(meta), env(basecall_model), env(kinetics)  , emit: meta_ext
    tuple val(meta), env(decode_with_fasta)              , emit: decode_with_fasta
    tuple val(meta), env(aligned_input)                  , emit: aligned_input
    tuple val(meta), env(m6a_calls)                      , emit: m6a_calls
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ont = meta.platform == 'ont'
    def m6a_check = m6aCallsCheck(bam)
    basecall_model = ''
    kinetics = ''
    decode_with_fasta = ''
    aligned_input = ''
    m6a_calls = ''
    // aligned_input: 'true' when the input (BAM or CRAM) has @SQ lines. Such input is prepared before fibertools-rs,
    // which takes no reference to decode a CRAM and panics on the SEQ-less secondary alignments minimap2 writes.
    // decode_with_fasta: 'true' when an aligned CRAM shares no contig name with the reference under a different M5,
    // so tools that would otherwise not find its reference may be given the pipeline FASTA. Given a mismatched
    // --reference, an embedded-reference CRAM decodes to wrong bases with exit 0, so any conflict means 'false'.
    def cram_ref_check = fasta ? """
    if [ -s input_sq.txt ] && ! samtools dict ${fasta} | grep '^@SQ' | awk -F'\\t' '
        function tag(t,   i) { for (i = 2; i <= NF; i++) if (substr(\$i, 1, 3) == t ":") return substr(\$i, 4); return "" }
        FNR == NR { m5[tag("SN")] = tag("M5"); next }
        (tag("SN") in m5) && tag("M5") != "" && tag("M5") != m5[tag("SN")] { bad = 1 }
        END { exit !bad }' - input_sq.txt; then
        decode_with_fasta=true
    fi
    """ : ''
    """
    export basecall_model="${basecall_model}"
    export kinetics="${kinetics}"
    decode_with_fasta=false
    aligned_input=false
    samtools view -H ${bam} | grep '^@SQ' > input_sq.txt || true
    if [ -s input_sq.txt ]; then
        aligned_input=true
    fi
    ${cram_ref_check}
    ${m6a_check}
    if [ $ont = 'true' ]; then
        basecall_model=\$(samtools view -H "${bam}" ${args} | awk -F'basecall_model=' '/basecall_model=/ {print \$2; exit}' | awk '{print \$1}' | tr -d '[:space:]')
    else
        kinetics=\$(samtools view -H ${bam} | awk '/--keep-kinetics/ {found=1} END {print (found ? "true" : "false")}')
        basecall_model="hifi_revio"
    fi
    """

    stub:
    // A real model name, so the Clair3/ClairS model lookups resolve and -stub runs reach the callers
    def stub_model = meta.platform == 'ont' ? 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0' : 'hifi_revio'
    def ont = meta.platform == 'ont'
    def m6a_check = m6aCallsCheck(bam)
    // The real kinetics header check, so a -stub run can give a sample's BAMs different kinetics
    // (tests/kinetics_pairing.nf.test); an empty placeholder reads as 'false'
    """
    export basecall_model="${stub_model}"
    export kinetics=""
    if [ $ont != 'true' ]; then
        kinetics=\$( (samtools view -H ${bam} 2> /dev/null || true) | awk '/--keep-kinetics/ {found=1} END {print (found ? "true" : "false")}')
    fi
    decode_with_fasta=false
    # The real header check, so -stub runs route aligned input as real runs do; an empty placeholder reads as unaligned
    aligned_input=false
    samtools view -H ${bam} 2> /dev/null | grep '^@SQ' > input_sq.txt || true
    if [ -s input_sq.txt ]; then
        aligned_input=true
    fi
    ${m6a_check}
    """
}

// m6a_calls: 'true' when an MM tag of the first 1000 records has an m6A entry (A+a, T-a or ChEBI 28871), 'false'
// when none has, 'unknown' when the records could not all be read. Only aux fields are decoded, so no CRAM needs a reference.
def m6aCallsCheck(bam) {
    """
    m6a_calls=\$( (samtools head -h 0 -n 1000 --input-fmt-option required_fields=0x800 ${bam} || echo '#read_failed') | awk -F'\\t' '
        \$0 == "#read_failed" { failed = 1; next }
        { n++; for (i = 12; i <= NF; i++) if (\$i ~ /^M[Mm]:Z:/) {
            k = split(substr(\$i, 6), e, ";")
            for (j = 1; j <= k; j++) if (e[j] ~ /^(A[+]|T-)([a-z]*a[a-z]*|28871)[.?]?(,|\$)/) found = 1
        } }
        END { print (found ? "true" : (failed || !n ? "unknown" : "false")) }')
    """
}
