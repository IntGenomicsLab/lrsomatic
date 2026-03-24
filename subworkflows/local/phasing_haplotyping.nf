// Import modules
include { LONGPHASE_PHASE as LONGPHASE_PHASE_GERMLINE       } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_PHASE as LONGPHASE_PHASE_SOMATIC        } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG                                } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { LONGPHASE_MODCALL as LONGPHASE_MODCALL_GERMLINE   } from '../../modules/local/longphase/modcall/main.nf'
include { LONGPHASE_MODCALL as LONGPHASE_MODCALL_SOMATIC    } from '../../modules/local/longphase/modcall/main.nf'
include { SAMTOOLS_INDEX                                    } from '../../modules/nf-core/samtools/index/main.nf'
include { BCFTOOLS_CONCAT                                   } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                     } from '../../modules/nf-core/bcftools/sort/main'


workflow PHASING_HAPLOTYPING {
    take:
    tumor_normal_bams // [meta, bam, bai]
    germline_vcf
    somatic_vcf
    fasta
    fai

    main:

    // SPLIT INTO PAIRED AND TUMOR ONLY
    tumor_normal_bams
        .branch { meta, _bams, _bai ->
            paired:      meta.paired_data
            tumor_only: !meta.paired_data
        }
        .set { branched_bams }

    branched_bams.paired
        .set{ paired_ch }

    branched_bams.tumor_only
        .map { meta, bam, bai ->
                def new_meta = meta.subMap('id',
                                        'paired_data',
                                        'platform',
                                        'sex',
                                        'fiber',
                                        'clair3_model',
                                        'clairS_model',
                                        'clairSTO_model',
                                        'kinetics')
                    return [ new_meta, bam, bai ]
            }
        .set{ tumor_only_ch }

    paired_ch
        .branch { meta, _bam, _bai ->
            normal: meta.type == "normal"
            tumor:  meta.type == "tumor"
        }
        .set {paired_ch_branched}

    paired_ch_branched.normal
        .map { meta, bam, bai ->
                def new_meta = meta.subMap('id',
                                        'paired_data',
                                        'platform',
                                        'sex',
                                        'fiber',
                                        'clair3_model',
                                        'clairS_model',
                                        'clairSTO_model',
                                        'kinetics')
                    return [ new_meta, bam, bai ]
            }
        .set{ paired_normal_ch }

    paired_ch_branched.tumor
        .map { meta, bam, bai ->
                def new_meta = meta.subMap('id',
                                        'paired_data',
                                        'platform',
                                        'sex',
                                        'fiber',
                                        'clair3_model',
                                        'clairS_model',
                                        'clairSTO_model',
                                        'kinetics')
                    return [ new_meta, bam, bai ]
            }
        .set{ paired_tumor_ch }

    tumor_only_ch
        .mix(paired_normal_ch)
        .set { normal_bams_w_tumoronly_ch }
    tumor_only_ch
        .mix(paired_tumor_ch)
        .set{ tumor_bams_ch}

    // MODCALL

    if (!params.skip_modcall) {

        LONGPHASE_MODCALL_GERMLINE (
            normal_bams_w_tumoronly_ch,
            fasta,
            fai
        )

        LONGPHASE_MODCALL_SOMATIC (
            tumor_bams_ch,
            fasta,
            fai
        )

    }
    // PHASING
    if (!params.skip_modcall) {
        normal_bams_w_tumoronly_ch
            .join(germline_vcf)
            .join(LONGPHASE_MODCALL_GERMLINE.out.mod_vcf)
            .map { meta, bam, bai, vcf, _tbi, mods->
                def svs = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_germline_input_ch }

        germline_vcf
            .join(somatic_vcf)
            .map { meta, germline_vcf, germline_tbi, somatic_vcf, somatic_tbi ->
                    def vcfs = [somatic_vcf, germline_vcf]
                    def tbis = [somatic_tbi, germline_tbi]
                    return [ meta, vcfs, tbis]
            }
            .set{germline_somatic_vcfs}
        BCFTOOLS_CONCAT(germline_somatic_vcfs)
        BCFTOOLS_CONCAT.out.vcf
                .set{concat_out}
        BCFTOOLS_SORT(concat_out)
        BCFTOOLS_SORT.out.vcf
            .set{germline_somatic_vcfs}

        tumor_bams_ch
            .join(germline_somatic_vcfs)
            .join(LONGPHASE_MODCALL_SOMATIC.out.mod_vcf)
            .map { meta, bam, bai, vcf, mods->
                def svs = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_somatic_input_ch }
    }
    else {
        normal_bams_w_tumoronly_ch
            .join(germline_vcf)
            .map { meta, bam, bai, vcf, _tbi ->
                def svs = []
                def mods = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_germline_input_ch }

        tumor_bams_ch
            .join(germline_somatic_vcfs)
            .join(LONGPHASE_MODCALL_SOMATIC.out.mod_vcf)
            .map { meta, bam, bai, vcf ->
                def svs = []
                def mods = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_somatic_input_ch }
    }

    LONGPHASE_PHASE_GERMLINE (
        longphase_phase_germline_input_ch,
        fasta,
        fai
    )

    LONGPHASE_PHASE_GERMLINE.out.snv_vcf
        .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf_index)
        .set{ phased_germline_vcf }

    LONGPHASE_PHASE_SOMATIC (
        longphase_phase_somatic_input_ch,
        fasta,
        fai
    )

    LONGPHASE_PHASE_SOMATIC.out.snv_vcf
        .join(LONGPHASE_PHASE_SOMATIC.out.snv_vcf_index)
        .set{ phased_germline_vcf }

    // HAPLOTAGING
    // remove type for merging


    if(!params.skip_modcall) {

        LONGPHASE_MODCALL_GERMLINE.out.mod_vcf
            .map { meta, mods ->
                def new_meta = meta.subMap('id',
                                    'paired_data',
                                    'platform',
                                    'sex',
                                    'fiber',
                                    'clair3_model',
                                    'clairS_model',
                                    'clairSTO_model',
                                    'kinetics')
                return [ new_meta, mods ]
            }
            .set{modcall_vcf_ch}

            tumor_only_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ tumor_only_ch }

            paired_tumor_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_tumor_ch }

            paired_normal_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "normal"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_normal_ch }

    }
    else {

            tumor_only_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ tumor_only_ch }

            paired_tumor_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_tumor_ch }

            paired_normal_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "normal"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_normal_ch }

    }

    tumor_only_ch
        .mix(paired_tumor_ch)
        .mix(paired_normal_ch)
        .set {longphase_haplotag_input_ch}

    LONGPHASE_HAPLOTAG (
        longphase_haplotag_input_ch,
        fasta,
        fai
    )

    LONGPHASE_HAPLOTAG.out.bam
        .set{ tumor_normal_hapbams_ch }

    SAMTOOLS_INDEX (
        tumor_normal_hapbams_ch
    )
    tumor_normal_hapbams_ch
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ tumor_normal_hapbams_ch }


    emit:
    tumor_normal_hapbams_ch
    phased_germline_vcf
}
