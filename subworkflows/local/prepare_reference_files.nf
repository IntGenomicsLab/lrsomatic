//
// Prepare reference files (unzipping and adding index)
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UNZIP as UNZIP_ALLELES         } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC              } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI            } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT              } from '../../modules/nf-core/unzip/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta
        ascat_alleles
        ascat_loci
        ascat_loci_gc
        ascat_loci_rt

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
        
        //
        // Prepare ASCAT files
        //
        
        // prepare ascat and controlfreec reference files
        if (!ascat_alleles) allele_files = Channel.empty()
        else if (ascat_alleles.endsWith(".zip")) {
            UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
            allele_files = UNZIP_ALLELES.out.unzipped_archive.map{ it[1] }
            ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
        } else allele_files = Channel.fromPath(ascat_alleles).collect()

        if (!ascat_loci) loci_files = Channel.empty()
        else if (ascat_loci.endsWith(".zip")) {
            UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
            loci_files = UNZIP_LOCI.out.unzipped_archive.map{ it[1] }
            ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
        } else loci_files = Channel.fromPath(ascat_loci).collect()

        if (!ascat_loci_gc) gc_file = Channel.value([])
        else if (ascat_loci_gc.endsWith(".zip")) {
            UNZIP_GC(Channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
            gc_file = UNZIP_GC.out.unzipped_archive.map{ it[1] }
            ch_versions = ch_versions.mix(UNZIP_GC.out.versions)
        } else gc_file = Channel.fromPath(ascat_loci_gc).collect()

        if (!ascat_loci_rt) rt_file = Channel.value([])
        else if (ascat_loci_rt.endsWith(".zip")) {
            UNZIP_RT(Channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
            rt_file = UNZIP_RT.out.unzipped_archive.map{ it[1] }
            ch_versions = ch_versions.mix(UNZIP_RT.out.versions)
        } else rt_file = Channel.fromPath(ascat_loci_rt).collect()

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai
        
        allele_files
        loci_files
        gc_file
        rt_file
        
        versions = ch_versions
}