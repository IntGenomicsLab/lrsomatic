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
        ch_versions = channel.empty()
        ch_prepared_fasta = channel.empty()
        allele_files = channel.empty()
        loci_files = channel.empty()
        gc_file = channel.empty()
        rt_file = channel.empty()

        // Check if fasta and gtf are zipped
        if (fasta.endsWith('.gz')){
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.file
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = [ [:], fasta ]
        }

        // if clair3 model is specified, then download that
        // otherwise use info in bam header and download that

        basecall_meta.map { meta, basecall_model_meta, kinetics_meta ->
            def id_new = basecall_model_meta ? clair3_modelMap.get(basecall_model_meta) : basecall_model_meta
            def meta_new = [id: id_new]
            def model = (!meta.clair3_model || meta.clair3_model.toString().trim() in ['', '[]']) ? clair3_modelMap.get(basecall_model_meta) : meta.clair3_model
            def download_prefix = ( basecall_model_meta == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
            def url = "${download_prefix}/${model}.tar.gz"
            return [ meta_new, url ]
        }
        .unique()
        .set{ clair3_model_urls }

        //
        // MODULE: Download model
        //

        WGET ( clair3_model_urls )

        ch_versions = ch_versions.mix(WGET.out.versions)

        //
        // MODULE: Untar model
        //

        UNTAR (
            WGET.out.outfile
        )

        UNTAR.out.untar.set { downloaded_clair3_models }

        //
        // MODULE: Index the fasta
        //

        SAMTOOLS_FAIDX (
            ch_prepared_fasta,
            [ [:], [] ],
            false
        )

        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai

        //
        // Prepare ASCAT files
        //

        // prepare ascat and controlfreec reference files
        if ( !params.skip_ascat ) {
            if (!ascat_alleles) allele_files = channel.empty()
            else if (ascat_alleles.endsWith(".zip")) {
                UNZIP_ALLELES(channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                allele_files = UNZIP_ALLELES.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
            } else allele_files = channel.fromPath(ascat_alleles).collect()

            if (!ascat_loci) loci_files = channel.empty()
            else if (ascat_loci.endsWith(".zip")) {
                UNZIP_LOCI(channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                loci_files = UNZIP_LOCI.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
            } else loci_files = channel.fromPath(ascat_loci).collect()

            if (!ascat_loci_gc) gc_file = channel.value([])
            else if ( ascat_loci_gc.endsWith(".zip") ) {
                UNZIP_GC(channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                gc_file = UNZIP_GC.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_GC.out.versions)
            } else gc_file = channel.fromPath(ascat_loci_gc).collect()

            if (!ascat_loci_rt) rt_file = channel.value([])
            else if (ascat_loci_rt.endsWith(".zip")) {
                UNZIP_RT(channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                rt_file = UNZIP_RT.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_RT.out.versions)
            } else rt_file = channel.fromPath(ascat_loci_rt).collect()
        }

    emit:
        prepped_fasta = ch_prepared_fasta
        prepped_fai = ch_prepared_fai

        allele_files
        loci_files
        gc_file
        rt_file
        downloaded_clair3_models

        versions = ch_versions
}
