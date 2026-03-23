// Import modules
include { LONGPHASE_PHASE           } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG        } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { LONGPHASE_MODCALL         } from '../../modules/local/longphase/modcall/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'
workflow PHASING_HAPLOTYPING {
    take:
    tumor_normal_bams // [meta, bam, bai]
    germline_vcf
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
        .set { longphase_modcall_input_ch }

    // MODCALL

    if (!params.skip_modcall) {

        LONGPHASE_MODCALL (
            longphase_modcall_input_ch,
            fasta,
            fai
        )

    }
    // PHASING
    if (!params.skip_modcall) {
        longphase_modcall_input_ch
            .join(germline_vcf)
            .join(LONGPHASE_MODCALL.out.mod_vcf)
            .map { meta, bam, bai, vcf, _tbi, mods->
                def svs = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_input_ch }
    }   
    else {
        longphase_modcall_input_ch
            .join(germline_vcf)
            .map { meta, bam, bai, vcf, _tbi->
                def svs = []
                def mods = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_input_ch }
    }
    LONGPHASE_PHASE (
        longphase_phase_input_ch,
        fasta,
        fai
    )

    LONGPHASE_PHASE.out.snv_vcf
        .join(LONGPHASE_PHASE.out.snv_vcf_index)
        .set{ phased_germline_vcf }

    // HAPLOTAGING
    // remove type for merging

    
    if(!params.skip_modcall) {

        LONGPHASE_MODCALL.out.mod_vcf
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
                .join(LONGPHASE_PHASE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods] 
                }
                .set{ tumor_only_ch }
            
            paired_tumor_ch
                .join(LONGPHASE_PHASE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods] 
                }
                .set{ paired_tumor_ch }
            
            paired_normal_ch
                .join(LONGPHASE_PHASE.out.snv_vcf)
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
                .join(LONGPHASE_PHASE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods] 
                }
                .set{ tumor_only_ch }
            
            paired_tumor_ch
                .join(LONGPHASE_PHASE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods] 
                }
                .set{ paired_tumor_ch }
            
            paired_normal_ch
                .join(LONGPHASE_PHASE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "normal"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods] 
                }
                .set{ paired_normal_ch }

    }

    tumor_only_ch
        .join(paired_tumor_ch)
        .join(paired_normal_ch)
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