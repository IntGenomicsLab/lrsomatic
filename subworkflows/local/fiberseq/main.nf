//
// Fiber-seq: m6A prediction (PacBio with kinetics), nucleosome and FIRE calling on the input reads (aligned input
// that fibertools-rs reads is cut to primary records as BAM first)
//

include { FIBERTOOLSRS_PREDICTM6A  } from '../../../modules/local/fibertoolsrs/predictm6a'
include { FIBERTOOLSRS_FIRE        } from '../../../modules/local/fibertoolsrs/fire'
include { FIBERTOOLSRS_NUCLEOSOMES } from '../../../modules/local/fibertoolsrs/nucleosomes'
include { FIBERTOOLSRS_QC          } from '../../../modules/local/fibertoolsrs/qc'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FIBERPREP } from '../../../modules/nf-core/samtools/view/main'
include { fibertoolsReadsInput     } from '../utils_nfcore_lrsomatic_pipeline'

workflow FIBERSEQ {

    take:
    ch_input_reads    // [meta_with_replicate, bam_or_cram] -- input reads, one per replicate (unaligned or aligned)
    ch_aligned_input  // [meta, 'true'|'false']  -- METAEXTRACT: the input has @SQ lines
    ch_decode_fasta   // [meta, 'true'|'false']  -- METAEXTRACT: decode this CRAM with --fasta (its M5s match)
    ch_fasta          // [[:], fasta]
    ch_fai            // [[:], fai]

    main:
    //
    // MODULE: SAMTOOLS_VIEW_FIBERPREP (nf-core samtools/view; 8 cpus, see modules.config) -- only for aligned input
    // that fibertools-rs reads.
    // fibertools-rs takes no reference, so it cannot decode a CRAM unless htslib finds one itself, and
    // `ft predict-m6a` panics on secondary alignments, which minimap2 writes with SEQ '*' by default. This keeps the
    // primary records (-F 0x900) as BAM, decoding a CRAM with --fasta when METAEXTRACT's M5 check allows it, and
    // from its embedded reference, UR: path or REF_PATH otherwise. MINIMAP2_ALIGN's samtools reset drops secondary
    // and supplementary records anyway. Every other sample goes on as given.
    // Input:  [meta, bam_or_cram, []]  -- aligned input of samples fibertoolsReadsInput() selects
    //         [meta, fasta, fai] or [[:], [], []]  -- the reference only where decode_with_fasta is 'true'
    // Output: .bam -- [meta, bam]  -- primary records only
    //
    // ch_input_reads' meta has gained the model fields since METAEXTRACT, so join on the fields both keep
    ch_input_reads
        .map { meta, bams -> [ [ meta.id, meta.type, meta.replicate ], meta, bams ] }
        .join(
            ch_aligned_input
                .join(ch_decode_fasta, failOnMismatch: true)
                .map { meta, aligned, decode -> [ [ meta.id, meta.type, meta.replicate ], aligned, decode ] },
            failOnMismatch: true
        )
        .branch { _key, meta, _bams, aligned, _decode ->
            prep: aligned == 'true' && fibertoolsReadsInput(meta, meta.kinetics)
            as_given: true
        }
        .set { ch_fiber_input }

    ch_fiber_input.prep
        .combine(ch_fasta)
        .combine(ch_fai)
        .multiMap { _key, meta, bams, _aligned, decode, meta_fasta, fasta, _meta_fai, fai ->
            reads: [ meta, bams, [] ]
            fasta: decode == 'true' ? [ meta_fasta, fasta, fai ] : [ [:], [], [] ]
        }
        .set { ch_fiberprep_in }

    SAMTOOLS_VIEW_FIBERPREP (
        ch_fiberprep_in.reads,
        ch_fiberprep_in.fasta,
        [ [:], [] ],
        [ [:], [] ],
        ''
    )

    ch_samplesheet = ch_fiber_input.as_given
        .map { _key, meta, bams, _aligned, _decode -> [ meta, bams ] }
        .mix(SAMTOOLS_VIEW_FIBERPREP.out.bam)
    // ch_samplesheet: [meta_with_replicate, bam_or_cram]  -- fibertools-rs input is now unaligned or primary-only BAM

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
