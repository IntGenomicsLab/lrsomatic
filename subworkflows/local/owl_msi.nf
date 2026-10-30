// IMPORT MODULES
include { OWL_PROFILE } from '../../modules/local/owl/profile/main'
include { OWL_SCORE   } from '../../modules/local/owl/score/main'

//
// Microsatellite instability with owl on the haplotagged tumour BAMs. owl is tumour-only (normals are not
// profiled) and licensed for PacBio data, so ONT samples run only with allow_ont.
//
workflow OWL_MSI {

    take:
    hapbams    // [meta, bam, bai]  -- haplotagged (HP/PS) BAMs, tumour and normal
    markers    // path: marker BED (chrom, start, end, motif)
    allow_ont  // boolean

    main:
    ch_tumor = hapbams.filter { meta, _bam, _bai ->
        meta.type == 'tumor' && (meta.platform == 'pb' || allow_ont)
    }

    OWL_PROFILE(ch_tumor, markers)
    OWL_SCORE(OWL_PROFILE.out.profile)

    emit:
    profile        = OWL_PROFILE.out.profile      // [meta, owl-profile.txt]
    scores         = OWL_SCORE.out.scores         // [meta, owl-scores.txt]
    motif_counts   = OWL_SCORE.out.motif_counts   // [meta, owl-motif-counts.txt]
    phased_sites   = OWL_SCORE.out.phased_sites   // [meta, owl-phased-sites.txt]
    phased_summary = OWL_SCORE.out.phased_summary // [meta, owl-phased-summary.txt]
}
