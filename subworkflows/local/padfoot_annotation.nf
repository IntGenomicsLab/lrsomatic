// IMPORT MODULES
include { PADFOOT as PADFOOT_SEVERUS_WAKHAN } from '../../modules/local/padfoot/main'
include { PADFOOT as PADFOOT_SAVANA         } from '../../modules/local/padfoot/main'
include { WGET as PADFOOT_WGET              } from '../../modules/nf-core/wget/main'
include { UNTAR as PADFOOT_UNTAR            } from '../../modules/nf-core/untar/main'

//
// Padfoot annotation of somatic SVs + CNAs, once per caller pair that produced output for a sample:
// Severus SVs + the top-ranked Wakhan integer-CN VCF, and SAVANA SVs + SAVANA absolute CN.
// Pass channel.empty() for a caller that did not run. Padfoot is not on bioconda: its source tree
// comes from params.padfoot_url (GitHub archive) or a local checkout in params.padfoot_dir.
//
workflow PADFOOT_ANNOTATION {

    take:
    severus_vcf      // [meta, severus_somatic.vcf.gz]
    wakhan_vcf_files // [meta, [wakhan_cna_*.vcf, ...]]  -- every fitted solution
    savana_vcf       // [meta, classified.somatic.vcf]
    savana_cna       // [meta, segmented_absolute_copy_number.tsv]
    fasta            // [[:], fasta]
    fai              // [[:], fai]
    annot            // [[:], padfoot_genome, gff | [], rm | []]

    main:
    ch_versions = channel.empty()

    if (params.padfoot_dir) {
        padfoot_src = channel.value([[id: 'padfoot'], file(params.padfoot_dir, type: 'dir', checkIfExists: true)])
    }
    else {
        PADFOOT_WGET( channel.value([[id: 'padfoot'], params.padfoot_url]) )
        PADFOOT_UNTAR( PADFOOT_WGET.out.outfile )
        padfoot_src = PADFOOT_UNTAR.out.untar
        ch_versions = ch_versions.mix(PADFOOT_WGET.out.versions)
    }
    // padfoot_src: [meta, dir]  -- padfoot.py + beds/

    //
    // MODULE: PADFOOT_SEVERUS_WAKHAN (label: process_medium)
    // Input:  [meta, severus_somatic.vcf.gz, 'severus', wakhan_cna_integers.vcf, 'wakhan']
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
        .map { meta, sv, cna -> [meta, sv, 'severus', cna, 'wakhan'] }
        .set { severus_wakhan_input }

    PADFOOT_SEVERUS_WAKHAN( severus_wakhan_input, fasta, fai, padfoot_src, annot )
    ch_versions = ch_versions.mix(PADFOOT_SEVERUS_WAKHAN.out.versions)

    //
    // MODULE: PADFOOT_SAVANA (label: process_medium)
    // Input:  [meta, classified.somatic.vcf, 'savana', segmented_absolute_copy_number.tsv, 'savana']
    //         SAVANA CN is only present when a fit was found, so the join drops unfitted samples
    //
    savana_vcf
        .join(savana_cna)
        .map { meta, sv, cna -> [meta, sv, 'savana', cna, 'savana'] }
        .set { savana_input }

    PADFOOT_SAVANA( savana_input, fasta, fai, padfoot_src, annot )
    ch_versions = ch_versions.mix(PADFOOT_SAVANA.out.versions)

    emit:
    versions = ch_versions // [versions.yml]
}
