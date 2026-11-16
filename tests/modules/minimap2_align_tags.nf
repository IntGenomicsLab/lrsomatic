// Test helper: wrapped in a workflow so the '.*:MINIMAP2_ALIGN' selector of conf/modules.config applies
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'

workflow MINIMAP2_ALIGN_PIPELINE_CONFIG {
    take:
    ch_reads   // [meta, bam]: meta as the pipeline sets it (id, type, replicate, platform, aligned_input, basecall_model)
    ch_fasta   // [meta2, fasta]

    main:
    MINIMAP2_ALIGN(ch_reads, ch_fasta, 'bam', 'bai', false, false)

    emit:
    bam = MINIMAP2_ALIGN.out.bam
}

process TAGGED_BAM_FIXTURE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam), path(fasta)   // minimap2-aligned BAM (carries tp, cm, s1, ...) and its reference

    output:
    tuple val(meta), path("tagged.bam")  , emit: bam
    tuple val(meta), path("reference.fa"), emit: fasta

    script:
    """
    zcat -f ${fasta} > reference.fa
    # 200 primary records plus the phasing, kinetics and read tags of a published haplotagged PacBio file
    samtools view -h -F 0x900 ${bam} \\
        | awk -F'\\t' -v OFS='\\t' '/^@/ { print; next }
            n < 200 { print \$0, "HP:i:1", "PS:i:1000", "PQ:i:30", "fi:B:C,1,2,3", "fp:B:C,4,5,6", "ri:B:C,1,2,3",
                "rp:B:C,4,5,6", "fn:i:5", "rn:i:4", "rq:f:0.999"; n++ }' \\
        | samtools view -b -o tagged.bam -
    if [ "\$(samtools view -c tagged.bam)" -ne 200 ]; then
        echo "tagged.bam: expected 200 records" >&2
        exit 1
    fi
    """
}
