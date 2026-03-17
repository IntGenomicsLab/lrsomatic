include { CLAIR3                    } from '../../modules/local/clair3/main.nf'
include { LONGPHASE_PHASE           } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG        } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'
include { CLAIRS                    } from '../../modules/local/clairs/main.nf'
include { BCFTOOLS_CONCAT           } from '../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT             } from '../../modules/nf-core/bcftools/sort'

include { DEEPVARIANT                                     } from '../../subworkflows/nf-core/deepvariant/main.nf'
include { DEEPSOMATIC                                     } from '../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS as GERMLINE_CONSENSUS   } from '../../subworkflows/local/small_variant_consensus.nf'
include { SMALL_VARIANT_CONSENSUS as SOMATIC_CONSENSUS    } from '../../subworkflows/local/small_variant_consensus.nf'



workflow TUMOR_NORMAL_HAPPHASE {
    take:
    mixed_bams
    fasta
    fai
    downloaded_clair3_models

    main:

    ch_versions = channel.empty()
    tumor_normal_severus = channel.empty()
    somatic_vep = channel.empty()
    germline_vep = channel.empty()

    // Branch input bams in normal and tumour
    mixed_bams
        .branch{ meta, _bam, _bai ->
            normal: meta.type == "normal"
            tumor: meta.type == "tumor"
        }
        .set{ mixed_bams }

    // Get normal bams and add platform/model info for Clair3 usage
    // remove type from so that information can be merged easier later

    downloaded_clair3_models
        .map{ meta, file ->
            def clair3_model = meta.id
            return [meta, clair3_model, file]
        }
        .set{downloaded_clair3_models}

    mixed_bams.normal
        .map{ meta, bam, bai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            clair3_model: meta.clair3_model,
                            clairS_model: meta.clairS_model,
                            clairSTO_model: meta.clairSTO_model,
                            kinetics: meta.kinetics]
            return [ new_meta, meta.clair3_model, bam, bai ]
        }
        .set { normal_bams_model }
    // [meta, clair3_model_id, bam, bai]  -- keyed by model ID for .combine() with downloaded_clair3_models

    normal_bams_model
        .combine(downloaded_clair3_models,by:1)
        .map {_clair3_model, meta_bam, bam, bai, _meta_model, model ->
            def platform = (meta_bam.platform == 'pb') ? 'hifi' : meta_bam.platform
            return [meta_bam, bam, bai, model, platform]
        }
        .set{ normal_bams }
    // [meta, bam, bai, clair3_model_dir, platform]  -- type excluded from meta; platform is "hifi" for PacBio

    /*
    .map{ basecall_model, meta, bam, bai, meta2, model ->
            def platform = (meta.platform == "pb") ? "hifi" : "ont"
            return [meta, bam, bai, model, platform]
        }
    */

    // Get tumour bams
    // remove type from so that information can be merged easier later
    mixed_bams.tumor
        .map{ meta, bam, bai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            clair3_model: meta.clair3_model,
                            clairS_model: meta.clairS_model,
                            clairSTO_model: meta.clairSTO_model,
                            kinetics: meta.kinetics]
            return[new_meta, bam, bai]
        }
        .set{ tumor_bams }
    // [meta, bam, bai]  -- type excluded from meta for downstream groupTuple merge

    //
    // MODULE: CLAIR3
    // small germline variant calling

    CLAIR3 (
        normal_bams,
        fasta,
        fai
    )
    
    normal_bams
        .map {meta, bam, bai, _model, _platform ->
            def intervals = []
            return [meta, bam, bai, intervals]
        }
    .set{deepvar_normal_bams}

    DEEPVARIANT (
        deepvar_normal_bams,
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

    CLAIR3.out.vcf
        .join(CLAIR3.out.tbi)
        .map { meta, vcf , tbi ->
            def new_meta = meta + [caller:'clair3']
            return [new_meta, vcf, tbi]
        }
        .set{clair3_ch}
        // [meta,deepvar_vcf,deepvar_index,clair3_vcf,clair3_index]
    
    clair3_ch
        .mix(deepvariant_ch)
        .set{mixed_vcfs}
    
    GERMLINE_CONSENSUS(
        mixed_vcfs,
        fasta,
        fai,
        params.germline_var_keep
    )
    

    // Add germline vcf to normal bams
    // remove clair3 model information
    normal_bams
        .join(GERMLINE_CONSENSUS.out.vcf)
        .map { meta, bam, bai, _clair3_model, _platform, vcf ->
            def svs = []
            def mods = []
            return [meta, bam, bai, vcf, svs, mods]
        }
        .set{ normal_bams_germlinevcf }
    // [meta, bam, bai, germline_vcf, [], []]  -- svs and mods are empty placeholders for LONGPHASE_PHASE input

    GERMLINE_CONSENSUS.out.vcf
        .map { meta, vcf ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { germline_vep }
    // [meta, clair3_vcf, []]  -- germline small variants for VEP annotation

    //
    // MODULE: LONGPHASE_PHASE
    //
    // Phase normals

    LONGPHASE_PHASE (
        normal_bams_germlinevcf,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

    // Add phased vcf to normal bams
    // Add type information back
    // both are needed for mixing with the tumor bams

    normal_bams
        .join(LONGPHASE_PHASE.out.snv_vcf)
        .map { meta, bam, bai, _clair3_model, _platform, vcf ->
            def new_meta = meta + [type: "normal"]
            def svs = []
            def mods = []
            return[new_meta, bam, bai, vcf, svs, mods]
        }
        .set{ normal_bams }
    // [meta+{type:"normal"}, bam, bai, phased_vcf, [], []]  -- type re-added; svs and mods are empty placeholders for LONGPHASE_HAPLOTAG

    // Add phased vcf to tumour bams and type information
    // mix with the normal bams

    tumor_bams
        .join(LONGPHASE_PHASE.out.snv_vcf)
        .map { meta, bam, bai, vcf ->
            def new_meta = meta + [type: "tumor"]
            def svs = []
            def mods = []
            return [new_meta, bam, bai, vcf, svs, mods]
        }
        .mix(normal_bams)
        .set{ mixed_bams_vcf }
    // [meta+{type}, bam, bai, phased_normal_vcf, [], []]  -- tumor and normal items both carry the same phased normal VCF

    //
    // MODULE: LONGPHASE_HAPLOTAG
    //

    // haplotag tumor and normal bams with normal vcf files for both
    LONGPHASE_HAPLOTAG (
        mixed_bams_vcf,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    // Get final tagged bams
    LONGPHASE_HAPLOTAG.out.bam
        .set{ mixed_hapbams }
    // [meta+{type}, haplotagged_bam]

    //
    // MODULE: SAMTOOLS_INDEX
    //
    // index the haplotaged bams

    SAMTOOLS_INDEX (
        mixed_hapbams
    )

    // Add index to channel
    mixed_bams_vcf
        .join(mixed_hapbams)
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ mixed_hapbams }
    // [meta+{type}, orig_bam, orig_bai, vcf, svs, mods, hapbam, hapbai]

    // Group everything back together in one channel
    mixed_hapbams
        .map { meta, _bam, _bai, _vcf, _snvs, _mods, hapbam, hapbai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            clair3_model: meta.clair3_model,
                            clairS_model: meta.clairS_model,
                            clairSTO_model: meta.clairSTO_model,
                            kinetics: meta.kinetics]
            return[new_meta, [[type: meta.type], hapbam], [[type: meta.type], hapbai]]
        }
        .groupTuple(size: 2)
        .map{ meta, bam, bai ->
            def normal_bam = bam[0][0].type == "normal" ? bam[0][1] : bam[1][1]
            def tumor_bam = bam[0][0].type == "tumor" ? bam[0][1] : bam[1][1]
            def normal_bai = bai[0][0].type == "normal" ? bai[0][1] : bai[1][1]
            def tumor_bai = bai[0][0].type == "tumor" ? bai[0][1] : bai[1][1]
            // Return channel
            return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
        }
        .join(LONGPHASE_PHASE.out.snv_vcf)
        .join(LONGPHASE_PHASE.out.snv_vcf_index)
        .set{tumor_normal_severus}
    // [meta, tumor_hapbam, tumor_bai, normal_hapbam, normal_bai, phased_vcf, phased_tbi]

    // Get ClairS input channel
    tumor_normal_severus
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _vcf, _tbi ->
            return[meta , tumor_bam, tumor_bai, normal_bam, normal_bai, meta.clairS_model]
        }
        .set { clairs_input }

    tumor_normal_severus
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _vcf, _tbi ->
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

    //
    // MODULE: CLAIRS
    //

    CLAIRS (
        clairs_input,
        fasta,
        fai
    )

    CLAIRS.out.vcfs
        .join(CLAIRS.out.tbi)
        .set{clairs_out}

    //
    // MODULE: BCFTOOLS_CONCAT
    //

    BCFTOOLS_CONCAT (
        clairs_out
    )

    //
    // MODULE: BCFTOOLS_SORT
    //

    BCFTOOLS_SORT (
        BCFTOOLS_CONCAT.out.vcf
    )
    BCFTOOLS_SORT.out.vcf.view()
    BCFTOOLS_SORT.out.tbi.view()

    DEEPSOMATIC.out.vcf
        .join(DEEPSOMATIC.out.vcf_index)
        .map{ meta, vcf, tbi ->
            def new_meta = meta + [caller:'deepsomatic']
            return [new_meta, vcf, tbi]
        }
        .set{deepsomatic_ch}

    BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi)
        .map { meta, vcf , tbi ->
            def new_meta = meta + [caller:'clairs']
            return [new_meta, vcf, tbi]
        }
        .set{clairs_ch}
        // [meta,deepvar_vcf,deepvar_index,clair3_vcf,clair3_index]
    clairs_ch.view()
    clairs_ch
        .mix(deepsomatic_ch)
        .set{mixed_somatic_vcfs}
    mixed_somatic_vcfs.view()
    SOMATIC_CONSENSUS(
        mixed_somatic_vcfs,
        fasta,
        fai,
        params.somatic_var_keep
    )

    SOMATIC_CONSENSUS.out.vcf
        .map { meta, vcf ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { somatic_vep }
        
    // [meta, sorted_clairs_vcf, []]  -- somatic small variants (SNV+indel merged) for VEP annotation

    emit:
    tumor_normal_severus
    somatic_vep
    germline_vep
    versions = ch_versions

}
