// IMPORT MODULES
include { RECONPLOT as RECONPLOT_ASCAT_SEVERUS  } from '../../modules/local/reconplot/main'
include { RECONPLOT as RECONPLOT_WAKHAN_SEVERUS } from '../../modules/local/reconplot/main'
include { RECONPLOT as RECONPLOT_SAVANA         } from '../../modules/local/reconplot/main'
include { WGET as RECONPLOT_PKG_WGET            } from '../../modules/nf-core/wget/main'
include { UNTAR as RECONPLOT_PKG_UNTAR          } from '../../modules/nf-core/untar/main'

//
// ReConPlot rearrangement + copy-number figures for every CN/SV caller pair that produced output for
// a sample: ASCAT + Severus, Wakhan + Severus, and SAVANA on its own. Pass channel.empty() for a
// caller that did not run. The wrapper (assets/reconplot) is shipped with the pipeline; the ReConPlot R
// package is staged as source from params.reconplot_pkg_url or a local checkout in params.reconplot_pkg_dir.
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
    if (params.reconplot_pkg_dir) {
        reconplot_pkg = channel.value([[id: 'reconplot_pkg'], file(params.reconplot_pkg_dir, type: 'dir', checkIfExists: true)])
    }
    else {
        RECONPLOT_PKG_WGET( channel.value([[id: 'reconplot_pkg'], params.reconplot_pkg_url]) )
        RECONPLOT_PKG_UNTAR( RECONPLOT_PKG_WGET.out.outfile )
        reconplot_pkg = RECONPLOT_PKG_UNTAR.out.untar
        ch_versions = ch_versions.mix(RECONPLOT_PKG_WGET.out.versions)
    }
    // reconplot_src: [meta, dir] -- wrapper;  reconplot_pkg: [meta, dir] -- R package source (conda installs it at run time)

    severus_sv_files = severus_vcf.map { meta, vcf -> [meta, [vcf]] }
    // severus_sv_files: [meta, [severus_somatic.vcf.gz]]

    //
    // MODULE: RECONPLOT_ASCAT_SEVERUS (label: process_low)
    // Input:  [meta, 'ascat', [segments.txt, purityploidy.txt, *BAF.txt], 'severus', [vcf]]
    //         the wrapper picks <sample>.tumour_tumourBAF.txt from the BAF tables by name
    //
    ascat_segments
        .join(ascat_purityploidy)
        .join(ascat_bafs)
        .map { meta, seg, pp, bafs -> [meta, [seg, pp, bafs].flatten()] }
        .join(severus_sv_files)
        .map { meta, cn, sv -> [meta, 'ascat', cn, 'severus', sv] }
        .set { ascat_input }

    RECONPLOT_ASCAT_SEVERUS( ascat_input, reconplot_src, reconplot_pkg, genome )
    ch_versions = ch_versions.mix(RECONPLOT_ASCAT_SEVERUS.out.versions)

    //
    // MODULE: RECONPLOT_WAKHAN_SEVERUS (label: process_low)
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

    RECONPLOT_WAKHAN_SEVERUS( wakhan_input, reconplot_src, reconplot_pkg, genome )
    ch_versions = ch_versions.mix(RECONPLOT_WAKHAN_SEVERUS.out.versions)

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

    RECONPLOT_SAVANA( savana_input, reconplot_src, reconplot_pkg, genome )
    ch_versions = ch_versions.mix(RECONPLOT_SAVANA.out.versions)

    emit:
    versions = ch_versions // [versions.yml]
}
