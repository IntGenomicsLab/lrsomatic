//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta

    main:
        ch_versions = Channel.empty()

        // Check if fasta and gtf are zipped

        //
        ch_prepared_fasta = Channel.empty()
        if (fasta.endsWith('.gz')){
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.file
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = [ [:], fasta ]
        }

        //
        // MODULE: Index the fasta
        //
        SAMTOOLS_FAIDX( ch_prepared_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ])
        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        versions = ch_versions
}