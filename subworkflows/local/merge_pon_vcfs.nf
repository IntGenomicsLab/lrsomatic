include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT   } from '../../modules/nf-core/bcftools/sort/main'

workflow MERGE_PON_VCFS {
    take:
    pon_vcfs  // channel: [meta, [vcf1, vcf2, ...]]  (bgzipped; indexes not required)

    main:
    ch_versions = Channel.empty()

    // BCFTOOLS_CONCAT auto-runs `tabix` on any VCF that lacks a .tbi/.csi index.
    // Pass an empty tbi list so it indexes every input file on the fly.
    ch_concat_input = pon_vcfs.map { meta, vcfs -> [meta, vcfs, []] }

    BCFTOOLS_CONCAT(ch_concat_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions_bcftools)

    // concat-with-overlaps output is not coordinate-sorted; sort and write final .tbi index
    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions_bcftools)

    emit:
    vcf_tbi  = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi) // [meta, vcf, tbi]
    versions = ch_versions
}
