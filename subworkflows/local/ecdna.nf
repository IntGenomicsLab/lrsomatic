//
// ecDNA and focal amplification: CoRAL reconstruction, AmpliconClassifier classification
//

// IMPORT MODULES
include { ASCAT_TO_CORAL_BED } from '../../modules/local/ascattocoralbed/main'
include { CORAL_SEED         } from '../../modules/local/coral/seed/main'
include { CORAL_RECONSTRUCT  } from '../../modules/local/coral/reconstruct/main'
include { CORAL_CYCLE        } from '../../modules/local/coral/cycle/main'
include { CORAL_PLOT         } from '../../modules/local/coral/plot/main'
include { AMPLICONCLASSIFIER } from '../../modules/local/ampliconclassifier/main'

workflow ECDNA {

    take:
    tumor_bam    // [meta, bam, bai]                 -- tumour BAMs only
    ascat_cnvs   // [meta, cnvs_txt]                 -- ASCAT.out.cnvs
    fai          // [[:], fai]                       -- value channel, reused by every sample
    data_repo    // [[:], dir]                       -- value channel; empty with --skip_ampliconclassifier
    coral_ref    // val 'hg38' | 't2t'
    ac_ref       // val 'GRCh38' | 'CHM13'

    main:
    //
    // MODULE: ASCAT_TO_CORAL_BED (label: process_single)
    // ASCAT writes cnvs.txt as chr/startpos/endpos/nMajor/nMinor; CoRAL wants
    // headerless BED with total CN last, spelled like the BAM's contigs.
    // Input:  [meta, cnvs_txt], [[:], fai]
    // Output: .bed -- [meta, bed]
    //
    ASCAT_TO_CORAL_BED (
        ascat_cnvs,
        fai
    )

    //
    // MODULE: CORAL_SEED (label: process_low)
    // Input:  [meta, cn_seg_bed, bam, bai]
    // Output: .seeds -- [meta, bed]  -- amplified intervals above --gain
    //
    // multiMap, not two reads of one channel: a channel feeding both a process and an
    // operator is not allowed, and reconstruct needs the same cn_seg/bam that seed used.
    ASCAT_TO_CORAL_BED.out.bed
        .join(tumor_bam, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, cn_seg, bam, bai ->
            seed: [meta, cn_seg, bam, bai]
            reconstruct: [meta, cn_seg, bam, bai]
            plot: [meta, bam, bai]
        }
        .set { coral_inputs }
    // coral_inputs.seed / .reconstruct: [meta, cn_seg_bed, bam, bai]; .plot: [meta, bam, bai]

    CORAL_SEED (
        coral_inputs.seed,
        coral_ref
    )

    //
    // A sample with no amplification above --gain yields an empty seed BED, and
    // CoRAL reconstruct errors on one. Size, not countLines(), so the file is
    // not staged just to be measured.
    //
    CORAL_SEED.out.seeds
        .branch { _meta, bed ->
            seeded: bed.toFile().length() > 0
            unseeded: true
        }
        .set { branched_seeds }

    branched_seeds.unseeded
        .subscribe { meta, _bed -> log.info("No amplified intervals found for ${meta.id}: skipping ecDNA reconstruction.") }

    //
    // MODULE: CORAL_RECONSTRUCT (label: process_high)
    // Input:  [meta, seeds, cn_seg_bed, bam, bai]
    // Output: .reconstruction -- [meta, dir]  -- graph/cycles/summary, named as AC expects
    //
    // remainder: false and no failOnMismatch: an unseeded sample is dropped here by design
    branched_seeds.seeded
        .join(coral_inputs.reconstruct, failOnDuplicate: true)
        .set { coral_reconstruct_input }
    // coral_reconstruct_input: [meta, seeds, cn_seg_bed, bam, bai]

    CORAL_RECONSTRUCT (
        coral_reconstruct_input
    )

    //
    // MODULE: CORAL_CYCLE (label: process_high) -- opt-in cycle re-extraction
    // Input:  [meta, reconstruction_dir]
    // Output: .reconstruction -- [meta, dir]  -- re-extracted cycles beside the copied graphs
    //
    if (params.coral_run_cycle) {
        CORAL_CYCLE (
            CORAL_RECONSTRUCT.out.reconstruction
        )
        ch_for_classifier = CORAL_CYCLE.out.reconstruction
    }
    else {
        ch_for_classifier = CORAL_RECONSTRUCT.out.reconstruction
    }
    // ch_for_classifier: [meta, dir]

    //
    // MODULE: CORAL_PLOT (label: process_medium)
    // Runs beside the classifier rather than in front of it: plotting is cosmetic
    // and must never gate classification.
    //
    ch_plots = channel.empty()
    if (params.coral_plot) {
        CORAL_RECONSTRUCT.out.reconstruction
            .join(coral_inputs.plot, failOnDuplicate: true)
            .set { coral_plot_input }
        // coral_plot_input: [meta, reconstruction_dir, bam, bai]

        CORAL_PLOT (
            coral_plot_input,
            coral_ref
        )
        ch_plots = CORAL_PLOT.out.plots
    }

    //
    // MODULE: AMPLICONCLASSIFIER (label: process_medium)
    // Input:  [meta, reconstruction_dir], [[:], data_repo], ac_ref
    // Output: .classification -- [meta, tsv]  -- amplicon_classification_profiles.tsv
    //
    ch_classification = channel.empty()
    if (!params.skip_ampliconclassifier) {
        AMPLICONCLASSIFIER (
            ch_for_classifier,
            data_repo,
            ac_ref
        )
        ch_classification = AMPLICONCLASSIFIER.out.classification
    }

    emit:
    classification = ch_classification                    // [meta, tsv]
    reconstruction = CORAL_RECONSTRUCT.out.reconstruction // [meta, dir]
    plots          = ch_plots                             // [meta, [files]]
}
