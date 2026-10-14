//
// Prepare reference files (unzipping and adding index)
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UNZIP as UNZIP_ALLELES         } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC              } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI            } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT              } from '../../modules/nf-core/unzip/main'
include { CLAIR3_MODEL                    } from '../../modules/local/clair3_model/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta           // str: path to reference FASTA (may be .gz)
        ascat_alleles   // str: path to ASCAT allele files (directory or .zip), or null
        ascat_loci      // str: path to ASCAT loci files (directory or .zip), or null
        ascat_loci_gc   // str: path to ASCAT GC correction file (.zip or direct), or null
        ascat_loci_rt   // str: path to ASCAT RT correction file (.zip or direct), or null
        clairsto_cna    // bool: also unzip the loci/allele/GC set for ClairS-TO's Verdict module
        basecall_meta   // [meta, basecall_model_str, kinetics_str]  -- from METAEXTRACT per sample
        clair3_modelMap // Map<basecall_model_str, clair3_model_name>  -- header basecall model to Clair3 model

    main:
        ch_versions = channel.empty()
        ch_prepared_fasta = channel.empty()
        allele_files = channel.empty()
        loci_files = channel.empty()
        gc_file = channel.empty()
        rt_file = channel.empty()

        // Decompress FASTA if gzipped; pass through as-is if already uncompressed
        if (fasta.endsWith('.gz')){
            //
            // MODULE: UNZIP_FASTA (PIGZ_UNCOMPRESS alias; label: process_medium)
            // Input:  [[:], fasta.gz]
            // Output: .file -- [[:], fasta]  -- decompressed FASTA
            //
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.file
        } else {
            ch_prepared_fasta = channel.value([ [:], fasta ])
        }
        // ch_prepared_fasta: [[:], fasta_path]  -- empty meta; uncompressed FASTA

        // Clair3 models: explicit clair3_model beats the BAM-header model. Clair3 v2 only reads PyTorch models;
        // the image bundles HKU's models and ONT's recent R10.4.1 ones, anything else is downloaded from
        // HKU's PyTorch conversions (ONT's Rerio catalogue for r9/r10 names)
        def clair3_bundled = [
            'hifi', 'hifi_revio', 'hifi_sequel2', 'ilmn', 'ont', 'ont_guppy5',
            'r1041_e82_400bps_hac_v410', 'r1041_e82_400bps_hac_v500', 'r1041_e82_400bps_hac_v520',
            'r1041_e82_400bps_hac_v520_with_mv', 'r1041_e82_400bps_hac_v600', 'r1041_e82_400bps_hac_v600_with_mv',
            'r1041_e82_400bps_hac_with_mv', 'r1041_e82_400bps_sup_v410', 'r1041_e82_400bps_sup_v430_bacteria_finetuned',
            'r1041_e82_400bps_sup_v500', 'r1041_e82_400bps_sup_v520', 'r1041_e82_400bps_sup_v520_with_mv',
            'r1041_e82_400bps_sup_with_mv', 'r941_prom_hac_g360+g422', 'r941_prom_sup_g5014'
        ]
        basecall_meta
            .map { meta, basecall_model_meta, _kinetics_meta ->
                def model = (!meta.clair3_model || meta.clair3_model.toString().trim() in ['', '[]']) ? clair3_modelMap.get(basecall_model_meta) : meta.clair3_model
                // Key on the model itself; keying on the header name let an explicit clair3_model add a duplicate
                return model
            }
            .unique()  // one entry per Clair3 model needed across all samples
            .branch { model ->
                bundled: model in clair3_bundled
                    return [ [id: model], [] ]
                download: true
                    def collection = model ==~ /^r(9|10).*/ ? 'clair3_models_rerio_pytorch' : 'clair3_models_pytorch'
                    return [ [id: model], "https://www.bio8.cs.hku.hk/clair3/${collection}/${model}" ]
            }
            .set { clair3_model_sources }
        // clair3_model_sources.bundled:  [meta(id=clair3_model_name), []]       -- read from the image
        // clair3_model_sources.download: [meta(id=clair3_model_name), url_str]  -- model directory URL

        //
        // MODULE: CLAIR3_MODEL (label: process_single)
        // Input:  [meta, url_str]  -- model name (id) + model directory URL
        // Output: .model -- [meta, model_dir]  -- pileup.pt and full_alignment.pt
        //
        CLAIR3_MODEL ( clair3_model_sources.download )

        clair3_model_sources.bundled
            .mix(CLAIR3_MODEL.out.model)
            .set { clair3_models }
        // clair3_models: [meta(id=clair3_model_name), model_dir or []]  -- [] = bundled with the Clair3 image

        //
        // MODULE: SAMTOOLS_FAIDX (label: process_single)
        // Input:  [[:], fasta, []]  -- empty meta + empty regions file (index full FASTA)
        //         false             -- do not write fai to stdout
        // Output: .fai -- [[:], fai_path]
        //
        SAMTOOLS_FAIDX (
            ch_prepared_fasta.map { meta, fa -> [meta, fa, []] },
            false
        )

        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai
        // ch_prepared_fai: [[:], fai_path]  -- empty meta

        //
        // ASCAT references: .zip or plain path each, emitted as flat file lists.
        // Loci/allele/GC are shared with Verdict; the RT file is ASCAT's alone (Verdict runs GC-only).
        if ( !params.skip_ascat || clairsto_cna ) {
            // Allele files: per-chromosome SNP allele frequency files (used for LogR/BAF calculation)
            if (!ascat_alleles) allele_files = channel.empty()
            else if (ascat_alleles.endsWith(".zip")) {
                // MODULE: UNZIP_ALLELES (UNZIP alias; label: process_single)
                // Input:  [meta(id=basename), [zip_file]]  -- collected zip
                // Output: .unzipped_archive -- [meta, dir]  -- flatMap lists the files inside
                UNZIP_ALLELES(channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                allele_files = UNZIP_ALLELES.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // allele_files: [path, path, ...]  -- all per-chromosome allele files collected
            } else allele_files = channel.fromPath(ascat_alleles).collect()

            // Loci files: per-chromosome SNP loci positions
            if (!ascat_loci) loci_files = channel.empty()
            else if (ascat_loci.endsWith(".zip")) {
                // MODULE: UNZIP_LOCI (UNZIP alias; label: process_single)
                UNZIP_LOCI(channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                loci_files = UNZIP_LOCI.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // loci_files: [path, path, ...]  -- all per-chromosome loci files collected
            } else loci_files = channel.fromPath(ascat_loci).collect()

            // GC correction file: genome-wide GC content per locus (optional)
            if (!ascat_loci_gc) gc_file = channel.value([])
            else if ( ascat_loci_gc.endsWith(".zip") ) {
                // MODULE: UNZIP_GC (UNZIP alias; label: process_single)
                UNZIP_GC(channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                gc_file = UNZIP_GC.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // gc_file: [path, ...]  -- GC correction file(s) collected
            } else gc_file = channel.fromPath(ascat_loci_gc).collect()
        }

        if ( !params.skip_ascat ) {
            // Replication timing correction file: RT correction per locus (optional)
            if (!ascat_loci_rt) rt_file = channel.value([])
            else if (ascat_loci_rt.endsWith(".zip")) {
                // MODULE: UNZIP_RT (UNZIP alias; label: process_single)
                UNZIP_RT(channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                rt_file = UNZIP_RT.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // rt_file: [path, ...]  -- RT correction file(s) collected
            } else rt_file = channel.fromPath(ascat_loci_rt).collect()
        }

    emit:
        prepped_fasta = ch_prepared_fasta  // [[:], fasta_path]  -- uncompressed reference FASTA
        prepped_fai   = ch_prepared_fai    // [[:], fai_path]    -- samtools FAI index

        // ASCAT reference files -- one flat list of paths each, no meta
        allele_files  // [path, ...]  -- per-chromosome allele frequency files
        loci_files    // [path, ...]  -- per-chromosome loci position files
        gc_file       // [path, ...]  -- GC correction file ([] if not provided)
        rt_file       // [path, ...]  -- replication timing correction file ([] if not provided)

        clair3_models  // [meta(id=clair3_model_name), model_dir or []]  -- [] = bundled with the Clair3 image

        versions = ch_versions
}
