//
// Fiber-seq: m6A prediction (PacBio with kinetics), nucleosome and FIRE calling on the unaligned BAMs
//

include { FIBERTOOLSRS_PREDICTM6A  } from '../../modules/local/fibertoolsrs/predictm6a'
include { FIBERTOOLSRS_FIRE        } from '../../modules/local/fibertoolsrs/fire'
include { FIBERTOOLSRS_NUCLEOSOMES } from '../../modules/local/fibertoolsrs/nucleosomes'
include { FIBERTOOLSRS_QC          } from '../../modules/local/fibertoolsrs/qc'

workflow FIBERSEQ {

    take:
    ch_samplesheet // [meta_with_replicate, bam] -- unaligned BAMs, one per replicate

    main:
    // Fiber-seq runs per replicate on the unaligned BAMs; replicates are merged after alignment
    if (!params.skip_normalfiber){
        // Process all samples (including normals) for fiber-seq
        ubams = ch_samplesheet
    }
    else {
        // Skip fiber-seq processing for normal samples; set aside normals to re-join later
        ch_samplesheet
        .branch { meta, _bams ->
            normal: meta.type == "normal"
            tumor: meta.type == "tumor"
            }
        .set { ch_ubams_normal_branching }
        // ch_ubams_normal_branching.normal: [meta, bam]  -- normal samples (held out)
        // ch_ubams_normal_branching.tumor:  [meta, bam]  -- tumor samples only

        normal_bams = ch_ubams_normal_branching.normal
        ubams = ch_ubams_normal_branching.tumor
    }
        // Branch by sequencing platform: PacBio needs m6A prediction, ONT does not
        ubams
        .branch{ meta, _bams ->
            pacBio: meta.platform == "pb"
            ont: meta.platform == "ont"
        }
        .set{ch_ubams_pacbio_ont_branching}
    // ch_ubams_pacbio_ont_branching.pacBio: [meta, bam]  -- PacBio samples
    // ch_ubams_pacbio_ont_branching.ont:    [meta, bam]  -- ONT samples (skip m6A)

    pacbio_bams = ch_ubams_pacbio_ont_branching.pacBio
    // Branch PacBio samples: only those with kinetics tags can have m6A predicted
    pacbio_bams
        .branch{meta, _bams ->
            kinetics: meta.kinetics == "true"
            noKinetics: meta.kinetics == "false"
        }
        .set{pacbio_bams}
    // pacbio_bams.kinetics:   [meta, bam]  -- PacBio with kinetics (mm/ml tags); m6A predictable
    // pacbio_bams.noKinetics: [meta, bam]  -- PacBio without kinetics; skip PREDICTM6A

    if (!params.skip_m6a) {
        //
        // MODULE: FIBERTOOLSRS_PREDICTM6A (label: process_high)
        // Input:  [meta, bam]  -- PacBio BAM with kinetics tags
        // Output: .bam -- [meta, bam]  -- BAM with m6A (MM/ML) tags added
        //
        FIBERTOOLSRS_PREDICTM6A (
            pacbio_bams.kinetics
        )
        // Merge PacBio with and without kinetics: both now have (or skip) m6A tags
        pacbio_bams.noKinetics
            .mix(FIBERTOOLSRS_PREDICTM6A.out.bam)
            .set{predicted_bams}
    }
    else {
        pacbio_bams.noKinetics
            .mix(pacbio_bams.kinetics)
            .set{predicted_bams}
    }
    // predicted_bams: [meta, bam]  -- all PacBio samples (m6A tags present where applicable)

    // Re-merge ONT and PacBio before fiber-seq branching
    ch_ubams_pacbio_ont_branching.ont
        .mix(predicted_bams)
        .set{fiber_branch}
    // fiber_branch (pre-split): [meta, bam]  -- all samples (ONT + PacBio, with m6A if applicable)

    // Branch on fiber-seq flag: only fiber-seq samples get nucleosome/FIRE calling
    fiber_branch
        .branch{ meta, _bams ->
            fiber: meta.fiber == "y"
            nonFiber: meta.fiber == "n"
        }
        .set{fiber_branch}
    // fiber_branch.fiber:    [meta, bam]  -- fiber-seq samples → nucleosome + FIRE calling
    // fiber_branch.nonFiber: [meta, bam]  -- non-fiber samples → passed through unchanged

    //
    // MODULE: FIBERTOOLSRS_NUCLEOSOMES (label: process_high)
    // Input:  [meta, bam]  -- fiber-seq BAM (with m6A tags for PacBio)
    // Output: .bam -- [meta, bam]  -- BAM with nucleosome footprint tags added
    //

    FIBERTOOLSRS_NUCLEOSOMES (
        fiber_branch.fiber
    )

    //
    // MODULE: FIBERTOOLSRS_FIRE (label: process_high)
    // Input:  FIBERTOOLSRS_NUCLEOSOMES.out.bam -- [meta, bam]  -- BAM with nucleosome tags
    // Output: .bam -- [meta, bam]  -- BAM with FIRE (Fiber-seq Inferred Regulatory Elements) tags
    //

    FIBERTOOLSRS_FIRE (
        FIBERTOOLSRS_NUCLEOSOMES.out.bam
    )

    if (!params.skip_normalfiber){
        // Re-merge fiber and non-fiber samples after FIRE annotation
        fiber_branch.nonFiber
        .mix(FIBERTOOLSRS_FIRE.out.bam)
        .set{ch_ubams}
    }
    else {
        // Re-merge fiber, non-fiber, and held-out normal samples
        fiber_branch.nonFiber
        .mix(normal_bams)
        .mix(FIBERTOOLSRS_FIRE.out.bam)
        .set{ch_ubams}
    }
    // ch_ubams (updated): [meta_with_replicate, bam]  -- all samples per-replicate;
    //   fiber-seq samples now carry nucleosome + FIRE tags; m6A tags for PacBio fiber-seq

    if(!params.skip_qc) {
        //
        // MODULE: FIBERTOOLSRS_QC (label: process_medium)
        // Input:  FIBERTOOLSRS_FIRE.out.bam -- [meta, bam]  -- annotated fiber-seq BAM
        // Output: QC reports for fiber-seq signal (written to outdir)
        //

        FIBERTOOLSRS_QC (
            FIBERTOOLSRS_FIRE.out.bam
        )
    }

    emit:
    bams = ch_ubams // [meta_with_replicate, bam] -- every input BAM; fiber-seq samples carry m6A/nucleosome/FIRE tags
}
