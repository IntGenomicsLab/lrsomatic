//
// Prepare reference files (unzipping and adding index)
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UNZIP as UNZIP_ALLELES         } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC              } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI            } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT              } from '../../modules/nf-core/unzip/main'
include { UNTAR                          } from '../../modules/nf-core/untar/main'
include { WGET                           } from '../../modules/nf-core/wget/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta
        ascat_alleles
        ascat_loci
        ascat_loci_gc
        ascat_loci_rt
        basecall_meta
        clair3_modelMap

    main:
        ch_versions = Channel.empty()
        ch_prepared_fasta = Channel.empty()
        allele_files = Channel.empty()
        loci_files = Channel.empty()
        gc_file = Channel.empty()
        rt_file = Channel.empty()

        // Check if fasta and gtf are zipped
        if (fasta.endsWith('.gz')){
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.file
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = [ [:], fasta ]
        }

        basecall_meta
            .map { meta, basecall_model_meta, kinetics_meta ->
            def meta_new = [id: basecall_model_meta]
            def model = clair3_modelMap.get(basecall_model_meta.toString().trim())
            def download_prefix = ( basecall_model_meta == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
            def url = "${download_prefix}/${model}.tar.gz"
            return [ meta_new, url ]
        }
        .unique()
        .set{ model_urls }

        //
        // MODULE: Download model
        //

        WGET ( model_urls )

        ch_versions = ch_versions.mix(WGET.out.versions)

        //
        // MODULE: Untar model
        //

        UNTAR (
            WGET.out.outfile
        )

        ch_versions = ch_versions.mix(UNTAR.out.versions)

        UNTAR.out.untar.set { downloaded_model_files }

        //
        // MODULE: Index the fasta
        //

        SAMTOOLS_FAIDX (
            ch_prepared_fasta,
            [ [:], "$projectDir/assets/dummy_file.txt" ],
            false
        )

        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

        //
        // Prepare ASCAT files
        //

        // prepare ascat and controlfreec reference files
        if ( !params.skip_ascat ) {
            if (!ascat_alleles) allele_files = Channel.empty()
            else if (ascat_alleles.endsWith(".zip")) {
                UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                allele_files = UNZIP_ALLELES.out.unzipped_archive.flatMap { it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
            } else allele_files = Channel.fromPath(ascat_alleles).collect()

            if (!ascat_loci) loci_files = Channel.empty()
            else if (ascat_loci.endsWith(".zip")) {
                UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                loci_files = UNZIP_LOCI.out.unzipped_archive.flatMap { it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
            } else loci_files = Channel.fromPath(ascat_loci).collect()

            if (!ascat_loci_gc) gc_file = Channel.value([])
            else if ( ascat_loci_gc.endsWith(".zip") ) {
                UNZIP_GC(Channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                gc_file = UNZIP_GC.out.unzipped_archive.flatMap { it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_GC.out.versions)
            } else gc_file = Channel.fromPath(ascat_loci_gc).collect()

            if (!ascat_loci_rt) rt_file = Channel.value([])
            else if (ascat_loci_rt.endsWith(".zip")) {
                UNZIP_RT(Channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                rt_file = UNZIP_RT.out.unzipped_archive.flatMap { it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_RT.out.versions)
            } else rt_file = Channel.fromPath(ascat_loci_rt).collect()
        }

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai

        allele_files
        loci_files
        gc_file
        rt_file
        downloaded_model_files

        versions = ch_versions
}
