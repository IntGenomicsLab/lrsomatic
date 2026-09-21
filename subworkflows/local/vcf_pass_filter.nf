// IMPORT MODULES
// Aliased so the software-versions report keeps this separate from the local
// BCFTOOLS_VIEW used by PHASING_HAPLOTYPING, which is a different bcftools build.
include { BCFTOOLS_VIEW as PASS_FILTER } from '../../modules/nf-core/bcftools/view/main'

//
// SUBWORKFLOW: VCF_PASS_FILTER
// Restrict a per-caller VCF to its PASS records before it is handed to the caller
// consensus, phasing, VEP and the report. DeepVariant and DeepSomatic emit a record for
// every site they evaluate (RefCall/GERMLINE/PON), and Clair3/ClairS keep their LowQual
// and NonSomatic calls, so without this the union is "every site every caller looked at".
// The per-caller VCFs published under variants/<caller>/ are produced elsewhere and are
// never filtered, so no calls are lost from the results directory.
//
// Include once per caller under an alias so conf/modules.config can give each one its own
// prefix, e.g. include { VCF_PASS_FILTER as DEEPVARIANT_PASS_FILTER }.
//
workflow VCF_PASS_FILTER {

    take:
    vcfs   // [meta, vcf, tbi]

    main:
    if (params.smallvar_filter_pass) {
        // The module declares `emit: index, optional: true`, so the index only exists because
        // conf/modules.config puts --write-index=tbi in ext.args for every alias. failOnMismatch
        // turns a missing index into an immediate error instead of silently dropping the sample.
        PASS_FILTER ( vcfs, [], [], [] )

        PASS_FILTER.out.vcf
            .join(PASS_FILTER.out.index, failOnMismatch: true, failOnDuplicate: true)
            .set{ filtered }
    }
    else {
        vcfs
            .set{ filtered }
    }

    emit:
    vcf = filtered   // [meta, vcf, tbi]
}
