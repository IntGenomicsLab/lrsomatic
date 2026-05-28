// IMPORT MODULES
include { CLAIR3                    } from '../../../modules/local/clair3/main.nf'

// IMPORT SUBWORKFLOWS
include { DEEPVARIANT                                     } from '../../../subworkflows/nf-core/deepvariant/main.nf'
include { SMALL_VARIANT_CONSENSUS as GERMLINE_CONSENSUS   } from '../../../subworkflows/local/small_variant_consensus.nf'

workflow PAIRED_SMALLVAR_GERMLINE {

    take:
    normal_bams   // [meta, normal_bam, normal_bai]  -- normal sample BAMs from T/N pairs
    fasta         // [[:], fasta]
    fai           // [[:], fai]
    clair3_models // [meta(id=model_name), model_dir]  -- downloaded Clair3 model directories

    main:
    germline_vcf = channel.empty()
    def germline_var_keep = params.germline_var_keep instanceof List ? params.germline_var_keep : params.germline_var_keep.toString().tokenize(',').collect { it.trim() }
    clair3_ch = channel.empty()
    deepvariant_ch = channel.empty()

    // COMBINE NORMAL BAMS WITH DOWNLOADED CLAIR3 MODELS
    // Clair3 requires the model directory path; models are keyed by model name (meta.id)
    if(germline_var_keep.contains('clair')) {

        // Extract model name from meta.id for combine-by key
        clair3_models
            .map{ meta, file ->
                def clair3_model_name = meta.id
                return [meta, clair3_model_name, file]
            }
            .set{clair3_models}
        // clair3_models: [meta(id=model_name), model_name_str, model_dir]

        // Emit [meta, clair3_model_name, bam, bai] to use model_name as the combine key
        normal_bams
            .map{ meta, bam, bai ->
                def new_meta = meta.subMap('id',
                                'paired_data',
                                'platform',
                                'sex',
                                'fiber',
                                'clair3_model',
                                'clairS_model',
                                'clairSTO_model',
                                'kinetics')
                return [ new_meta, meta.clair3_model, bam, bai ]
            }
            .set { normal_bams_model }
        // normal_bams_model: [meta, clair3_model_name, bam, bai]
        //   clair3_model_name is the join key used by .combine(clair3_models, by:1)

        //
        // MODULE: CLAIR3 (label: process_high)
        // Input:  [meta, bam, bai, model_dir, platform_str]
        //         fasta / fai
        // Output: .vcf -- [meta, vcf]  -- germline SNVs/indels
        //         .tbi -- [meta, tbi]
        //
        normal_bams_model
            .combine(clair3_models,by:1)  // join on clair3_model_name
            .map {_clair3_model, meta_bam, bam, bai, _meta_model, model ->
                def platform = (meta_bam.platform == 'pb') ? 'hifi' : meta_bam.platform
                return [meta_bam, bam, bai, model, platform]
            }
            .set{ clair3_input_ch }
        // clair3_input_ch: [meta, bam, bai, model_dir, platform_str]
        //   platform_str: 'hifi' for PacBio ('pb' → 'hifi'), otherwise meta.platform (e.g. 'ont')

        CLAIR3 (
            clair3_input_ch,
            fasta,
            fai
        )

        CLAIR3.out.vcf
            .join(CLAIR3.out.tbi)
            .map { meta, vcf , tbi ->
                def new_meta = meta + [caller:'clair3']
                return [new_meta, vcf, tbi]
            }
            .set{clair3_ch}
        // clair3_ch: [meta(+caller:'clair3'), vcf, tbi]
    }

    // DEEPVARIANT
    if(germline_var_keep.contains('deepvariant')) {

        //
        // SUBWORKFLOW: DEEPVARIANT (nf-core)
        // Input:  [meta, bam, bai, []]  -- [] is empty intervals (genome-wide)
        //         fasta / fai
        //         [[:],[]] x2  -- empty PAR/GFF interval files (not used for WGS)
        // Output: .vcf       -- [meta, vcf]
        //         .vcf_index -- [meta, tbi]
        //
        normal_bams
            .map {meta, bam, bai  ->
                def new_meta = meta.subMap('id',
                                    'paired_data',
                                    'platform',
                                    'sex',
                                    'fiber',
                                    'clair3_model',
                                    'clairS_model',
                                    'clairSTO_model',
                                    'kinetics')
                def intervals = []
                return [new_meta, bam, bai, intervals]
            }
        .set{deepvariant_input_ch}
        // deepvariant_input_ch: [meta, bam, bai, []]

        DEEPVARIANT (
            deepvariant_input_ch,
            fasta,
            fai,
            [[:],[]],  // PAR regions (not used)
            [[:],[]]   // GFF annotation (not used)
        )

        DEEPVARIANT.out.vcf
            .join(DEEPVARIANT.out.vcf_index)
            .map{ meta, vcf, tbi ->
                def new_meta = meta + [caller:'deepvariant']
                return [new_meta, vcf, tbi]
            }
            .set{deepvariant_ch}
        // deepvariant_ch: [meta(+caller:'deepvariant'), vcf, tbi]
    }

    // COMBINE GERMLINE VARIATION
    // If both callers requested: run consensus subworkflow; otherwise pass through single-caller output
    if (germline_var_keep.size() > 1) {
        // Mix both caller VCFs into a single channel for GERMLINE_CONSENSUS
        clair3_ch
            .mix(deepvariant_ch)
            .set{combined_germline_ch}
        // combined_germline_ch: [meta(+caller), vcf, tbi]  -- one item per caller per sample

        // SUBWORKFLOW: GERMLINE_CONSENSUS (SMALL_VARIANT_CONSENSUS alias)
        // Normalise, annotate with caller ID, intersect, and combine per params
        GERMLINE_CONSENSUS(
            combined_germline_ch,
            fasta,
            fai,
            params.prioritize_caller_germline,
            params.germline_var_combine
        )
        GERMLINE_CONSENSUS.out.vcf
            .join(GERMLINE_CONSENSUS.out.tbi)
            .set{ germline_vcf }
        // germline_vcf: [meta(+caller from consensus), vcf, tbi]
    }
    else if (germline_var_keep == ['clair']) {
        clair3_ch
            .set{germline_vcf}
    }
    else if (germline_var_keep == ['deepvariant']) {
        deepvariant_ch
            .set{germline_vcf}
    }

    // Strip 'caller' field from final germline VCF meta (not needed downstream)
    germline_vcf
        .map{ meta, vcf, tbi ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, vcf, tbi]
        }
        .set{germline_vcf}

    emit:
    germline_vcf  // [meta, vcf, tbi]  -- final germline VCF (Clair3, DeepVariant, or consensus)
}
