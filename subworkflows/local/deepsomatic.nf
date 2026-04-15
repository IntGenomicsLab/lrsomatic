include { DEEPSOMATIC_MAKEEXAMPLES        } from '../../modules/local/deepsomatic/makeexamples/main'
include { DEEPSOMATIC_CALLVARIANTS        } from '../../modules/local/deepsomatic/callvariants/main'
include { DEEPSOMATIC_POSTPROCESSVARIANTS } from '../../modules/local/deepsomatic/postprocessvariants/main'

workflow DEEPSOMATIC {
    take:
    ch_input     // [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
    //              normal_bam/bai may be [] for tumor-only mode
    _ch_intervals // [[:], []]  -- empty intervals (genome-wide calling)
    ch_fasta     // [[:], fasta]
    ch_fai       // [[:], fai]
    ch_gzi       // [[:], gzi]  -- bgzipped FASTA index (empty if FASTA is not bgzipped)

    main:

    //
    // MODULE: DEEPSOMATIC_MAKEEXAMPLES (label: process_high)
    // Input:  [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
    // Output: .examples -- [meta, [tfrecord shards...]]  -- serialised pileup examples
    //         .gvcf     -- [meta, [gvcf tfrecord shards...]]
    //
    DEEPSOMATIC_MAKEEXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi)

    //
    // MODULE: DEEPSOMATIC_CALLVARIANTS (label: process_gpu / process_high)
    // Input:  DEEPSOMATIC_MAKEEXAMPLES.out.examples -- [meta, [tfrecord shards...]]
    // Output: .call_variants_tfrecords -- [meta, tfrecord]  -- DNN variant call records
    //
    DEEPSOMATIC_CALLVARIANTS(DEEPSOMATIC_MAKEEXAMPLES.out.examples)

    // Join CALLVARIANTS output with MAKEEXAMPLES gVCF records only when generate_gvcf is true.
    // ch_postproc_input: [meta, call_tfrecord, [gvcf_tfrecords...], [], []]
    //   trailing [] are for optional candidate positions and haplotype outputs (unused)
    ch_postproc_input = params.generate_gvcf
        ? DEEPSOMATIC_CALLVARIANTS.out.call_variants_tfrecords.join(
              DEEPSOMATIC_MAKEEXAMPLES.out.gvcf, failOnMismatch: true
          ).map { meta, call_tfrecord, gvcf_tfrecords ->
              [meta, call_tfrecord, gvcf_tfrecords, [], []]
          }
        : DEEPSOMATIC_CALLVARIANTS.out.call_variants_tfrecords.map { meta, call_tfrecord ->
              [meta, call_tfrecord, [], [], []]
          }

    //
    // MODULE: DEEPSOMATIC_POSTPROCESSVARIANTS (label: process_medium)
    // Input:  [meta, call_tfrecord, [gvcf_tfrecords...], [], []]
    // Output: .vcf       -- [meta, vcf]   -- somatic variant calls (VCF)
    //         .vcf_index -- [meta, tbi]
    //         .gvcf      -- [meta, gvcf]  -- genome VCF (all sites)
    //         .gvcf_index-- [meta, tbi]
    //
    DEEPSOMATIC_POSTPROCESSVARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )

    emit:
    vcf        = DEEPSOMATIC_POSTPROCESSVARIANTS.out.vcf        // [meta, vcf]
    vcf_index  = DEEPSOMATIC_POSTPROCESSVARIANTS.out.vcf_index  // [meta, tbi]
    gvcf       = DEEPSOMATIC_POSTPROCESSVARIANTS.out.gvcf       // [meta, gvcf]
    gvcf_index = DEEPSOMATIC_POSTPROCESSVARIANTS.out.gvcf_index // [meta, tbi]
}
