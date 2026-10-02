// IMPORT MODULES
include { PADFOOT as PADFOOT_SEVERUS_WAKHAN } from '../../modules/local/padfoot/main'
include { PADFOOT as PADFOOT_SAVANA         } from '../../modules/local/padfoot/main'

//
// Padfoot annotation of somatic SVs + CNAs, once per caller pair that produced output for a sample:
// Severus SVs + the top-ranked Wakhan integer-CN VCF, and SAVANA SVs + SAVANA absolute CN, each with the
// caller's fitted purity/ploidy table so gene copy number is labelled against the tumour ploidy.
// Pass channel.empty() for a caller that did not run. Padfoot is not on bioconda: it ships inside the
// module's container (containers/padfoot/), so nothing is downloaded at run time.
//
workflow PADFOOT_ANNOTATION {

    take:
    severus_vcf            // [meta, severus_somatic.vcf.gz]
    wakhan_vcf_files       // [meta, [wakhan_cna_*.vcf, ...]]  -- every fitted solution
    wakhan_solutions_ranks // [meta, solutions_ranks.tsv]  -- purity/ploidy per solution (rank 1 = solution_1)
    savana_vcf             // [meta, classified.somatic.vcf]
    savana_cna             // [meta, segmented_absolute_copy_number.tsv]
    savana_purity_ploidy   // [meta, fitted_purity_ploidy.tsv]
    fasta                  // [[:], fasta]
    fai              // [[:], fai]
    annot            // [[:], padfoot_genome, gff | [], rm | []]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: PADFOOT_SEVERUS_WAKHAN (label: process_medium)
    // Input:  [meta, severus_somatic.vcf.gz, 'severus', wakhan_cna_integers.vcf, 'wakhan', solutions_ranks.tsv | []]
    //
    // Wakhan writes every fitted solution; solution_1/ holds the top-ranked one
    wakhan_vcf_files
        .map { meta, vcfs ->
            def integers = [vcfs].flatten().findAll { vcf -> vcf.name.endsWith('_wakhan_cna_integers.vcf') }
            return [meta, integers.find { vcf -> vcf.toString().contains('/solution_1/') } ?: integers[0]]
        }
        .filter { _meta, vcf -> vcf != null }
        .set { wakhan_best_cna }
    // wakhan_best_cna: [meta, wakhan_cna_integers.vcf]

    severus_vcf
        .join(wakhan_best_cna)
        .join(wakhan_solutions_ranks, remainder: true)       // the fit table is optional for Padfoot
        .filter { _meta, sv, _cna, _ranks -> sv != null }    // drops right-only rows of the remainder join
        .map { meta, sv, cna, ranks -> [meta, sv, 'severus', cna, 'wakhan', ranks ?: []] }
        .set { severus_wakhan_input }

    PADFOOT_SEVERUS_WAKHAN( severus_wakhan_input, fasta, fai, annot )
    ch_versions = ch_versions.mix(PADFOOT_SEVERUS_WAKHAN.out.versions)

    //
    // MODULE: PADFOOT_SAVANA (label: process_medium)
    // Input:  [meta, classified.somatic.vcf, 'savana', segmented_absolute_copy_number.tsv, 'savana', fitted_purity_ploidy.tsv | []]
    //         SAVANA CN is only present when a fit was found, so the join drops unfitted samples
    //
    savana_vcf
        .join(savana_cna)
        .join(savana_purity_ploidy, remainder: true)
        .filter { _meta, sv, _cna, _pp -> sv != null }
        .map { meta, sv, cna, pp -> [meta, sv, 'savana', cna, 'savana', pp ?: []] }
        .set { savana_input }

    PADFOOT_SAVANA( savana_input, fasta, fai, annot )
    ch_versions = ch_versions.mix(PADFOOT_SAVANA.out.versions)

    emit:
    versions = ch_versions // [versions.yml]
}
