// Test helper: wrapped in a workflow so the '.*:PREPARE_ALIGNED' selector of conf/modules.config applies
include { PREPARE_ALIGNED } from '../../modules/local/prepare_aligned/main'

workflow PREPARE_ALIGNED_PIPELINE_CONFIG {
    take:
    ch_reads       // [meta, bam]: meta as the pipeline sets it (id, type, replicate, platform, basecall_model)
    ch_fasta_fai   // [meta2, fasta, fai], a value channel

    main:
    PREPARE_ALIGNED(ch_reads, ch_fasta_fai, params.aligned_format)

    emit:
    alignment = PREPARE_ALIGNED.out.alignment
    index     = PREPARE_ALIGNED.out.index
}

process ALIGNED_FIXTURE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam), path(fasta)   // coordinate-sorted BAM and a reference (aligned to it or not)

    output:
    tuple val(meta), path("tagged.bam")                            , emit: bam
    tuple val(meta), path("reference.fa"), path("reference.fa.fai"), emit: fasta_fai

    script:
    """
    zcat -f ${fasta} > reference.fa
    samtools faidx reference.fa
    # The first 200 records, still sorted, with an old read group and the phasing and kinetics tags of a
    # published haplotagged PacBio file
    samtools view -h ${bam} \\
        | awk -F'\\t' -v OFS='\\t' '/^@RG/ { next } /^@/ { print; next }
            n == 0 { print "@RG", "ID:old", "SM:old" }
            n < 200 { sub(/\\tRG:Z:[^\\t]*/, ""); print \$0, "RG:Z:old", "HP:i:1", "PS:i:1000", "PQ:i:30", "PC:i:5",
                "fi:B:C,1,2,3", "fp:B:C,4,5,6", "ri:B:C,1,2,3", "rp:B:C,4,5,6", "fn:i:5", "rq:f:0.999"; n++ }' \\
        | samtools view -b -o tagged.bam -
    if [ "\$(samtools view -c tagged.bam)" -ne 200 ]; then
        echo "tagged.bam: expected 200 records" >&2
        exit 1
    fi
    """
}
