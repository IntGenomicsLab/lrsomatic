// IMPORT MODULES
include { CLAIR3                    } from '../../../modules/local/clair3/main.nf'

// IMPORT SUBWORKFLOWS
include { DEEPVARIANT                                     } from '../../../subworkflows/nf-core/deepvariant/main.nf'
include { SMALL_VARIANT_CONSENSUS as GERMLINE_CONSENSUS   } from '../../../subworkflows/local/small_variant_consensus.nf'

workflow PAIRED_SMALLVAR_GERMLINE {

    take:
    normal_bams // [ meta, normal_bam, normal_bai ]
    fasta
    fai
    clair3_models

    main:
    germline_vcf = channel.empty()
    germline_tbi = channel.empty()
    // COMBINE NORMAL BAMS WITH DOWNLOADED CLAIR3 MODELS
    if(params.germline_var_keep != 'deepvariant') {

        clair3_models
            .map{ meta, file ->
                def clair3_model_name = meta.id
                return [meta, clair3_model_name, file]
            }
            .set{clair3_models}
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

        // CLAIR3
        normal_bams_model
            .combine(clair3_models,by:1)
            .map {_clair3_model, meta_bam, bam, bai, _meta_model, model ->
                def platform = (meta_bam.platform == 'pb') ? 'hifi' : meta_bam.platform
                return [meta_bam, bam, bai, model, platform]
            }
            .set{ clair3_input_ch }

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
    }
    // DEEPVARIANT
    if(params.germline_var_keep != 'clair') {

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

        DEEPVARIANT (
            deepvariant_input_ch,
            fasta,
            fai,
            [[:],[]],
            [[:],[]]
        )

        DEEPVARIANT.out.vcf
            .join(DEEPVARIANT.out.vcf_index)
            .map{ meta, vcf, tbi ->
                def new_meta = meta + [caller:'deepvariant']
                return [new_meta, vcf, tbi]
            }
            .set{deepvariant_ch}
    }
    // COMBINE GERMLINE VARIATION
    if (params.germline_var_keep != 'clair' && params.germline_var_keep != 'deepvariant' ) {
        clair3_ch
            .mix(deepvariant_ch)
            .set{combined_germline_ch}

        GERMLINE_CONSENSUS(
            combined_germline_ch,
            fasta,
            fai,
            params.germline_var_keep
        )
        GERMLINE_CONSENSUS.out.vcf
            .join(GERMLINE_CONSENSUS.out.tbi)
            .set{ germline_vcf }
    }
    else if (params.germline_var_keep == 'clair') {
        clair3_ch
            .set{germline_vcf}
    }
    else if (params.germline_var_keep == 'deepvariant') {
        deepvariant_ch
            .set{germline_vcf}
    }

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
    germline_vcf
}
