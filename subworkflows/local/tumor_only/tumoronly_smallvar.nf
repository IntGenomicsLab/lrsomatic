// IMPORT MODULES
include { CLAIRSTO                  } from '../../../modules/local/clairsto/main.nf'
include { VCFSPLIT                  } from '../../../modules/local/vcfsplit/main.nf'

// IMPORT SUBWORKFLOWS
include { DEEPVARIANT                                   } from '../../../subworkflows/nf-core/deepvariant/main.nf'
include { DEEPSOMATIC                                   } from '../../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS as GERMLINE_CONSENSUS } from '../../../subworkflows/local/small_variant_consensus.nf'
include { SMALL_VARIANT_CONSENSUS as SOMATIC_CONSENSUS  } from '../../../subworkflows/local/small_variant_consensus.nf'


workflow TUMORONLY_SMALLVAR {

    take:
    tumor_bams // [ meta, tumor_bams, tumor_bai ]
    fasta
    fai
    pon_channel

    main:

    // empty channel emission

    ch_versions = channel.empty()
    somatic_vcf = channel.empty()
    germline_vcf = channel.empty()
    somatic_tbi = channel.empty()
    germline_tbi = channel.empty()

    // CLAIRS-TO (SOMATIC/NONGERMLINE VARIANT CALLING)

    if(params.somatic_var_keep != 'deepvariant') {
        tumor_bams
            .map { meta, bam, bai ->
                return [ meta, bam, bai, meta.clairSTO_model]
            }
            .combine(pon_channel)
            .set{ clairsto_input_ch}
        CLAIRSTO (
            clairsto_input_ch,
            fasta,
            fai
        )


        // SPLIT CLAIRSTO GERMLINE AND SOMATIC VARIATION

        CLAIRSTO.out.indel_vcf
                    .join(CLAIRSTO.out.snv_vcf)
                    .set{ clairsto_combined_vcf }

        VCFSPLIT (
            clairsto_combined_vcf
        )

        VCFSPLIT.out.germline_vcf
            .join(VCFSPLIT.out.germline_tbi)
            .map { meta, vcf, tbi ->
                def new_meta = meta + [caller:'clairs-to']
                return [ new_meta, vcf, tbi]
            }
            .set{clairsto_germline_ch}

        VCFSPLIT.out.somatic_vcf
            .join(VCFSPLIT.out.somatic_tbi)
            .map { meta, vcf, tbi ->
                def new_meta = meta + [caller:'clairs-to']
                return [ new_meta, vcf, tbi]
            }
            .set{clairsto_somatic_ch}
    }
    // DEEPVARIANT
    if(params.somatic_var_keep != 'clair') {
        tumor_bams
            .map { meta, bam, bai  ->
                def intervals = []
                return [meta,bam,bai, intervals]
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

    // COMBINE GERMLINE VARIANTS
    if (params.germline_var_keep != 'clair' | params.germline_var_keep != 'deepvariant' ) {
        clairsto_germline_ch
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
            .set{germline_vcf}
    }
    else if (params.germline_var_keep == 'clair') {
        clairsto_germline_ch
            .set{germline_vcf}
    }
    else if (params.germline_var_keep == 'deepvariant') {
        deepvariant_ch
            .set{germline_vcf}
    }
    // DEEPSOMATIC
    if(params.somatic_var_keep != 'clair') {
        tumor_bams
            .map { meta, tumor_bam, tumor_bai ->
                def normal_bam = []
                def normal_bai = []
                return [meta,normal_bam,normal_bai,tumor_bam,tumor_bai]
            }
            .set{deepsomatic_input_ch}

        DEEPSOMATIC (
            deepsomatic_input_ch,
            [[:],[]],
            fasta,
            fai,
            [[:],[]]
        )
        DEEPSOMATIC.out.vcf
            .join(DEEPSOMATIC.out.vcf_index)
            .map{ meta, vcf, tbi ->
                def new_meta = meta + [caller:'deepsomatic']
                return [new_meta, vcf, tbi]
            }
            .set{deepsomatic_ch}
    }
    // COMBINE SOMATIC VARIATION
    if (params.somatic_var_keep != 'clair' | params.somatic_var_keep != 'deepvariant' ) {
        clairsto_somatic_ch
            .mix(deepsomatic_ch)
            .set{combined_somatic_ch}

        SOMATIC_CONSENSUS(
            combined_somatic_ch,
            fasta,
            fai,
            params.somatic_var_keep
        )
        SOMATIC_CONSENSUS.out.vcf
            .join(SOMATIC_CONSENSUS.out.tbi)
            .set{somatic_vcf}
    }
    else if (params.somatic_var_keep == 'clair') {
        clairsto_somatic_ch
            .set{somatic_vcf}
    }
    else if (params.somatic_var_keep == 'deepvariant') {
        deepvariant_ch
            .set{somatic_vcf}
    }

    somatic_vcf
        .map{ meta, vcf, tbi  ->
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
        .set{somatic_vcf}

    germline_vcf
        .map{ meta, vcf, tbi  ->
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
    somatic_vcf
    germline_vcf


}
