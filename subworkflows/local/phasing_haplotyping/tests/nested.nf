// Test helper, not part of the pipeline. nf-test runs the workflow under test at the top level, so its processes
// would be named PHASING_HAPLOTYPING:<process> and the '.*:PHASING_HAPLOTYPING:...' selectors in conf/modules.config
// would not match them (REHEADER_GERMLINE and REJOIN_* would then write no index). One level of nesting, as in the
// pipeline, keeps the pipeline's process configuration.

include { PHASING_HAPLOTYPING } from '../main'

workflow NESTED_PHASING_HAPLOTYPING {
    take:
    tumor_normal_bams
    germline_vcf
    somatic_vcf
    fasta
    fai

    main:
    PHASING_HAPLOTYPING ( tumor_normal_bams, germline_vcf, somatic_vcf, fasta, fai )

    emit:
    tumor_normal_hapbams_ch = PHASING_HAPLOTYPING.out.tumor_normal_hapbams_ch
    phased_germline_vcf     = PHASING_HAPLOTYPING.out.phased_germline_vcf
    phased_somatic_vcf      = PHASING_HAPLOTYPING.out.phased_somatic_vcf
}
