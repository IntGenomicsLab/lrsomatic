// IMPORT MODULES
include { CLAIRS                    } from '../../../modules/local/clairs/main.nf'
include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort'

// IMPORT SUBWORKFLOWS
include { DEEPSOMATIC                                    } from '../../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS as SOMATIC_CONSENSUS   } from '../../../subworkflows/local/small_variant_consensus.nf'

workflow PAIRED_SMALLVAR_SOMATIC {

    take:
    tumor_normal_bams // [ meta, tumor_bam, tumor_bai, normal_hapbam, normal_bai ]
    fasta
    fai

    main:
    somatic_vcf = channel.empty()
    somatic_tbi = channel.empty()

    // CLAIRS
    if(params.somatic_var_keep != 'deepvariant') {
        tumor_normal_bams
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                return[meta , tumor_bam, tumor_bai, normal_bam, normal_bai, meta.clairS_model]
            }
            .set { clairs_input }

        CLAIRS (
            clairs_input,
            fasta,
            fai
        )

        // CONCAT CLAIRS INDEL AND SNV OUTPUT

        CLAIRS.out.vcfs
            .join(CLAIRS.out.tbi)
            .set{clairs_out}

        BCFTOOLS_CONCAT (
            clairs_out
        )

        BCFTOOLS_SORT (
            BCFTOOLS_CONCAT.out.vcf
        )

        BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT.out.tbi)
            .map { meta, vcf , tbi ->
                def new_meta = meta + [caller:'clairs']
                return [new_meta, vcf, tbi]
            }
            .set{clairs_ch}
    }
    // DEEPSOMATIC

    if(params.somatic_var_keep != 'clair') {

        tumor_normal_bams
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    return [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
            }
            .set{ deepsomatic_input }

        DEEPSOMATIC (
            deepsomatic_input,
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
    // COMBINE GERMLINE VARIATION
    if (params.somatic_var_keep != 'clair' && params.somatic_var_keep != 'deepvariant' ) {
        clairs_ch
            .mix(deepsomatic_ch)
            .set{combine_somatic_ch}

        SOMATIC_CONSENSUS(
            combine_somatic_ch,
            fasta,
            fai,
            params.somatic_var_keep
        )

        SOMATIC_CONSENSUS.out.vcf
            .join(SOMATIC_CONSENSUS.out.tbi)
            .set{ somatic_vcf }
    }
    else if (params.somatic_var_keep == 'clair') {
        clairs_ch
            .set{somatic_vcf}
    }
    else if (params.somatic_var_keep == 'deepvariant') {
        deepsomatic_ch
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
    emit:
    somatic_vcf
}
