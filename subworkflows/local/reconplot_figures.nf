// IMPORT MODULES
include { RECONPLOT as RECONPLOT_SEVERUS_ASCAT  } from '../../modules/local/reconplot/main'
include { RECONPLOT as RECONPLOT_SEVERUS_WAKHAN } from '../../modules/local/reconplot/main'
include { RECONPLOT as RECONPLOT_SAVANA         } from '../../modules/local/reconplot/main'

//
// ReConPlot rearrangement + copy-number figures for every CN/SV caller pair that produced output for
// a sample: ASCAT + Severus, Wakhan + Severus, and SAVANA on its own. Pass channel.empty() for a
// caller that did not run. The wrapper (assets/reconplot) ships with the pipeline; the ReConPlot R package
// ships inside the module's container (containers/reconplot/), so nothing is downloaded at run time.
//
workflow RECONPLOT_FIGURES {

    take:
    severus_vcf            // [meta, severus_somatic.vcf.gz]
    ascat_segments         // [meta, segments.txt]
    ascat_purityploidy     // [meta, purityploidy.txt]
    ascat_bafs             // [meta, [*BAF.txt]]
    wakhan_bed_files       // [meta, [bed_output/*.bed, ...]]  -- every fitted solution
    wakhan_solutions_ranks // [meta, solutions_ranks.tsv]
    savana_cna             // [meta, segmented_absolute_copy_number.tsv]
    savana_bedpe           // [meta, classified.somatic.bedpe]
    savana_purity_ploidy   // [meta, fitted_purity_ploidy.tsv]
    savana_allele_counts   // [meta, allele_counts_hetSNPs.bed]  -- absent without an SNP source
    genome                 // ReConPlot genome preset: hg38, hg19, T2T, mm10 or mm39

    main:
    ch_versions = channel.empty()

    reconplot_src = channel.value([[id: 'reconplot'], file("${projectDir}/assets/reconplot", type: 'dir', checkIfExists: true)])
    // reconplot_src: [meta, dir] -- the wrapper (run_reconplot.R + R/ + VERSION)

    severus_sv_files = severus_vcf.map { meta, vcf -> [meta, [vcf]] }
    // severus_sv_files: [meta, [severus_somatic.vcf.gz]]

    //
    // MODULE: RECONPLOT_SEVERUS_ASCAT (label: process_low)
    // Input:  [meta, 'ascat', [segments.txt, purityploidy.txt, *BAF.txt], 'severus', [vcf]]
    //         the wrapper picks <sample>.tumour_tumourBAF.txt from the BAF tables by name
    //
    // ASCAT writes an empty segments.txt (and NA purity/ploidy) when it finds no solution: nothing to draw for that pair
    ascat_segments
        .filter { _meta, seg -> seg.size() > 0 && seg.countLines() > 1 }
        .join(ascat_purityploidy)
        .join(ascat_bafs)
        .map { meta, seg, pp, bafs -> [meta, [seg, pp, bafs].flatten()] }
        .join(severus_sv_files)
        .map { meta, cn, sv -> [meta, 'ascat', cn, 'severus', sv] }
        .set { ascat_input }

    RECONPLOT_SEVERUS_ASCAT( ascat_input, reconplot_src, genome )
    ch_versions = ch_versions.mix(RECONPLOT_SEVERUS_ASCAT.out.versions)

    //
    // MODULE: RECONPLOT_SEVERUS_WAKHAN (label: process_low)
    // Input:  [meta, 'wakhan', [HP_1.bed, HP_2.bed, solutions_ranks.tsv], 'severus', [vcf]]
    //         the two allele-specific segment BEDs of the top-ranked solution (solution_1/)
    //
    wakhan_bed_files
        .map { meta, beds ->
            def hp = [beds].flatten().findAll { bed -> bed.name ==~ /.*_copynumbers_segments_HP_[12]\.bed/ }
            def best = hp.findAll { bed -> bed.toString().contains('/solution_1/') } ?: hp
            return [meta, best.unique { bed -> bed.name }]
        }
        .filter { _meta, beds -> beds.size() == 2 }
        .join(wakhan_solutions_ranks)
        .map { meta, beds, ranks -> [meta, beds + [ranks]] }
        .join(severus_sv_files)
        .map { meta, cn, sv -> [meta, 'wakhan', cn, 'severus', sv] }
        .set { wakhan_input }

    RECONPLOT_SEVERUS_WAKHAN( wakhan_input, reconplot_src, genome )
    ch_versions = ch_versions.mix(RECONPLOT_SEVERUS_WAKHAN.out.versions)

    //
    // MODULE: RECONPLOT_SAVANA (label: process_low)
    // Input:  [meta, 'savana', [cna.tsv, somatic.bedpe, fitted_purity_ploidy.tsv, allele_counts.bed], 'savana', []]
    //         single-source mode; allele_counts is optional, and a sample with allele counts but no CN fit
    //         only exists on the right of the remainder join ([meta, null, bed]) and is dropped
    //
    savana_cna
        .join(savana_bedpe)
        .join(savana_purity_ploidy)
        .join(savana_allele_counts, remainder: true)
        .filter { row -> row[1] != null }
        .map { meta, cna, bedpe, pp, hetsnp -> [meta, 'savana', [cna, bedpe, pp, hetsnp].findAll { f -> f != null }, 'savana', []] }
        .set { savana_input }

    RECONPLOT_SAVANA( savana_input, reconplot_src, genome )
    ch_versions = ch_versions.mix(RECONPLOT_SAVANA.out.versions)

    emit:
    versions = ch_versions // [versions.yml]
}
