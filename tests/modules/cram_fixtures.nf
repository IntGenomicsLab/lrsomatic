// Test helper, not part of the pipeline: builds the CRAM inputs of the CRAM-input module tests from an aligned test BAM
process CRAM_FIXTURES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam), path(fasta), path(other_fasta)   // BAM aligned to fasta; other_fasta shares its contig names

    output:
    tuple val(meta), path("aligned_embedded.cram"), emit: embedded
    tuple val(meta), path("aligned_noref.cram")   , emit: noref
    tuple val(meta), path("unaligned.cram")       , emit: unaligned
    tuple val(meta), path("reference.fa")         , emit: fasta
    tuple val(meta), path("other_reference.fa")   , emit: other_fasta
    // The BAM plus secondary alignments as minimap2 writes them (SEQ '*'), as BAM and as both CRAMs, and the BAM's
    // primary records alone: the records a filter dropping secondary and supplementary alignments must leave
    tuple val(meta), path("secondary.bam")         , emit: secondary_bam
    tuple val(meta), path("secondary_embedded.cram"), emit: secondary_embedded
    tuple val(meta), path("secondary_noref.cram")  , emit: secondary_noref
    tuple val(meta), path("primary_only.bam")      , emit: primary_only

    script:
    """
    zcat -f ${fasta} > reference.fa
    zcat -f ${other_fasta} > other_reference.fa
    # Secondary copies of the first 5 forward or reverse primary records: flag + 256, MAPQ 0, no SEQ/QUAL, no tags
    if [ "\$(samtools view -c -f 0x100 ${bam})" -ne 0 ]; then
        echo "${bam} already has secondary alignments; primary_only.bam below assumes it has none" >&2
        exit 1
    fi
    samtools view -h ${bam} \\
        | awk -F'\\t' -v OFS='\\t' '/^@/ { print; next } { print }
            n < 5 && (\$2 == 0 || \$2 == 16) { \$2 = \$2 + 256; \$5 = 0; \$10 = "*"; \$11 = "*"
                s = \$1; for (i = 2; i <= 11; i++) s = s OFS \$i; print s; n++ }' \\
        | samtools view -b -o secondary.bam -
    if [ "\$(samtools view -c -f 0x100 secondary.bam)" -ne 5 ]; then
        echo "secondary.bam: expected 5 injected secondary alignments" >&2
        exit 1
    fi
    samtools view -b -F 0x800 -o primary_only.bam ${bam}
    # The CRAMs are written against a copy that is then deleted, so their @SQ UR: tags point nowhere and
    # the *noref.cram files can only be decoded with a reference handed to samtools
    cp reference.fa ur_gone.fa
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 --output-fmt-option embed_ref=1 -o aligned_embedded.cram ${bam}
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 -o aligned_noref.cram ${bam}
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 --output-fmt-option embed_ref=1 -o secondary_embedded.cram secondary.bam
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 -o secondary_noref.cram secondary.bam
    rm ur_gone.fa ur_gone.fa.fai
    samtools reset -O cram,version=3.0 -o unaligned.cram ${bam}
    """
}
