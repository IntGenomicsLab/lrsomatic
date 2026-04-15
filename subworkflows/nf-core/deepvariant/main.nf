include { DEEPVARIANT_MAKEEXAMPLES        } from '../../../modules/nf-core/deepvariant/makeexamples/main'
include { DEEPVARIANT_MAKEEXAMPLES        } from '../../../modules/nf-core/deepvariant/makeexamples/main'
include { DEEPVARIANT_CALLVARIANTS        } from '../../../modules/nf-core/deepvariant/callvariants/main'
include { DEEPVARIANT_POSTPROCESSVARIANTS } from '../../../modules/nf-core/deepvariant/postprocessvariants/main'

workflow DEEPVARIANT {
    take:
    ch_input   // channel: [ val(meta), path(input), path(index), path(intervals)]
    ch_fasta   // channel: [ val(meta2), path(fasta) ]
    ch_fai     // channel: [ val(meta3), path(fai) ]
    ch_gzi     // channel: [ val(meta4), path(gzi) ]
    ch_par_bed // channel: [ val(meta5), path(par_bed) ]

    main:

    DEEPVARIANT_MAKEEXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi, ch_par_bed)

    DEEPVARIANT_CALLVARIANTS(DEEPVARIANT_MAKEEXAMPLES.out.examples)

    // Input to postprocessing step needs the variant calls from CALLVARIANTS,
    // and optionally the gvcfs from MAKEEXAMPLES (only when params.generate_gvcf is true).
    ch_intervals = ch_input.map { meta, _input, _index, intervals -> [ meta, intervals ] }

    ch_call_and_gvcf = params.generate_gvcf
        ? DEEPVARIANT_CALLVARIANTS.out.call_variants_tfrecords.join(
              DEEPVARIANT_MAKEEXAMPLES.out.gvcf, failOnMismatch: true
          )
        : DEEPVARIANT_CALLVARIANTS.out.call_variants_tfrecords.map { meta, tfrecord ->
              [meta, tfrecord, []]
          }

    ch_postproc_input = ch_call_and_gvcf.join(
        DEEPVARIANT_MAKEEXAMPLES.out.small_model_calls,
        failOnMismatch: true
    ).join(
        ch_intervals,
        failOnMismatch: true
    )

    DEEPVARIANT_POSTPROCESSVARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )

    emit:
    vcf        = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf
    vcf_index  = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf_index
    gvcf       = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf
    gvcf_index = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf_index
}
