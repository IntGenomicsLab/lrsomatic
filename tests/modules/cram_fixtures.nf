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

    script:
    """
    zcat -f ${fasta} > reference.fa
    zcat -f ${other_fasta} > other_reference.fa
    # The CRAMs are written against a copy that is then deleted, so their @SQ UR: tags point nowhere and
    # aligned_noref.cram can only be decoded with a reference handed to samtools
    cp reference.fa ur_gone.fa
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 --output-fmt-option embed_ref=1 -o aligned_embedded.cram ${bam}
    samtools view -C --reference ur_gone.fa --output-fmt-option version=3.0 -o aligned_noref.cram ${bam}
    rm ur_gone.fa ur_gone.fa.fai
    samtools reset -O cram,version=3.0 -o unaligned.cram ${bam}
    """
}
