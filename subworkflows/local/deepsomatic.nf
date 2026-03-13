include { DEEPSOMATIC_MAKEEXAMPLES        } from '../../modules/local/deepsomatic/makeexamples/main'
include { DEEPSOMATIC_CALLVARIANTS        } from '../../modules/local/deepsomatic/callvariants/main'
include { DEEPSOMATIC_POSTPROCESSVARIANTS } from '../../modules/local/deepsomatic/postprocessvariants/main'

workflow DEEPSOMATIC {
    take:
    ch_input   // channel: [ val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index)]
    ch_intervals
    ch_fasta   // channel: [ val(meta2), path(fasta) ]
    ch_fai     // channel: [ val(meta3), path(fai) ]
    ch_gzi     // channel: [ val(meta4), path(gzi) ]

    main:

    DEEPSOMATIC_MAKEEXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi)

    DEEPSOMATIC_CALLVARIANTS(DEEPSOMATIC_MAKEEXAMPLES.out.examples)

    // Input to postprocessing step needs both the gvcfs from MAKEEXAMPLES and the variant
    // calls from CALLVARIANTS. Joining on meta, which is assumed to be unique.
   

    ch_postproc_input = DEEPSOMATIC_CALLVARIANTS.out.call_variants_tfrecords.join(
        DEEPSOMATIC_MAKEEXAMPLES.out.gvcf,
        failOnMismatch: true
    ).map { meta, call_tfrecord, gvcf_tfrecords ->
        [meta, call_tfrecord, gvcf_tfrecords, [], []]
    }

    DEEPSOMATIC_POSTPROCESSVARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )

    emit:
    vcf        = DEEPSOMATIC_POSTPROCESSVARIANTS.out.vcf
    vcf_index  = DEEPSOMATIC_POSTPROCESSVARIANTS.out.vcf_index
    gvcf       = DEEPSOMATIC_POSTPROCESSVARIANTS.out.gvcf
    gvcf_index = DEEPSOMATIC_POSTPROCESSVARIANTS.out.gvcf_index
}
