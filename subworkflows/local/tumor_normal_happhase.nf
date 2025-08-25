include { CLAIR3 } from '../../modules/local/clair3/main.nf'
include { LONGPHASE_PHASE } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'
include { CLAIRS                    } from '../../modules/local/clairs/main.nf'
include {ENSEMBLVEP_VEP as SOMATIC_VEP} from '../../modules/nf-core/ensemblvep/vep/main.nf'
include {ENSEMBLVEP_VEP as GERMLINE_VEP} from '../../modules/nf-core/ensemblvep/vep/main.nf'

workflow TUMOR_NORMAL_HAPPHASE {
    take:
    mixed_bams
    fasta
    fai
    clair3_modelMap
    clairs_modelMap
    downloaded_model_files

    main:

    ch_versions = Channel.empty()

    // Branch input bams in normal and tumour
    mixed_bams
        .branch{ meta, bam, bai ->
            normal: meta.type == "normal"
            tumor: meta.type == "tumor"
        }
        .set{ mixed_bams }

    // Get normal bams and add platform/model info for Clair3 usage
    // remove type from so that information can be merged easier later

    downloaded_model_files
        .map{ meta, file ->
            def basecall_model = meta.id
            return [basecall_model, meta, file]
        }
        .set{downloaded_model_files}

     mixed_bams.normal
        .map{ meta, bam, bai ->
            def basecall_model = meta.basecall_model
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            basecall_model: meta.basecall_model]
            return [ basecall_model, new_meta, bam, bai ]
        }
        .set { normal_bams_model }

    normal_bams_model
        .combine(downloaded_model_files,by:0)
        .map{ basecall_model, meta, bam, bai, meta2, model ->
            def platform = (meta.platform == "pb") ? "hifi" : "ont"
            def clair3_model = (!meta.clair3_model || meta.clair3_model.toString().trim() in ['', '[]']) ? clair3_modelMap.get(meta.basecall_model) : meta.clair3_model
            return [meta, bam, bai, model, platform]
        }
        .set{ normal_bams }

    // normal_bams -> meta:         [id, paired_data, platform, sex, fiber, basecall_model]
    //                bam:          list of concatenated aligned bams
    //                bai:          indexes for bam files
    //                clair3_model: clair3 model name
    //                platform:     name of sequencing platform


    // Get tumour bams
    // remove type from so that information can be merged easier later
    mixed_bams.tumor
        .map{ meta, bam, bai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            basecall_model: meta.basecall_model]
            return[new_meta, bam, bai]
        }
        .set{ tumor_bams }

    // tumor_bams -> meta:  [id, paired_data, platform, sex, fiber, basecall_model]
    //                bam:  list of concatenated aligned bams
    //                bai:  indexes for bam files

    //
    // MODULE: CLAIR3
    //
    // small germline variant calling

    CLAIR3 (
        normal_bams,
        fasta,
        fai
    )

    

    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    // Add germline vcf to normal bams
    // remove clair3 model information

    normal_bams
        .join(CLAIR3.out.vcf)
        .map { meta, bam, bai, clair3_model, platform, vcf ->
            def svs = []
            def mods = []
            return [meta, bam, bai, vcf, svs, mods]
        }
        .set{ normal_bams_germlinevcf }

    // normal_bams -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                bam:  list of concatenated aligned bams
    //                bai:  indexes for bam files
    //                vcf:  normal small germline variant vcf
    //                svs:  structural variant vcf (empty)
    //                mods: modcall-generated VCF with modifications (empty)


    GERMLINE_VEP (
        CLAIR3.out.vcf,
        params.genome,
        "homo_sapiens",
        111,
        '',
        fasta,
        []
    )

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
        .join(LONGPHASE_PHASE.out.vcf)
        .map { meta, bam, bai, clair3_model, platform, vcf ->
            def new_meta = meta + [type: "normal"]
            def snvs = []
            def mods = []
            return[new_meta, bam, bai, vcf, snvs, mods]
        }
        .set{ normal_bams }

    // normal_bams -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                bam:  list of concatenated aligned bams
    //                bai:  indexes for bam files
    //                vcf:  normal small germline variant vcf
    //                svs:  structural variant vcf (empty)
    //                mods: modcall-generated VCF with modifications (empty)


    // Add phased vcf to tumour bams and type information
    // mix with the normal bams
    tumor_bams
        .join(LONGPHASE_PHASE.out.vcf)
        .map { meta, bam, bai, vcf ->
            def new_meta = meta + [type: "tumor"]
            def snvs = []
            def mods = []
            return [new_meta, bam, bai, vcf, snvs, mods]
        }
        .mix(normal_bams)
        .set{ mixed_bams_vcf }
    // mixed_bams_vcf -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                   bam:  list of concatenated aligned bams
    //                   bai:  indexes for bam files
    //                   vcf:  normal small germline variant vcf
    //                   svs:  structural variant vcf (empty)
    //                   mods: modcall-generated VCF with modifications (empty)

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

    // mixed_hapbams -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                  bams: haplotagged aligned bams


    //
    // MODULE: SAMTOOLS_INDEX
    //
    // index the haplotaged bams

    SAMTOOLS_INDEX (
        mixed_hapbams
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Add index to channel
    mixed_bams_vcf
        .join(mixed_hapbams)
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ mixed_hapbams }

    // mixed_hapbams -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                  bams: haplotagged aligned bams
    //                  bais: indexes for bam files

    // Group everything back together in one channel
    mixed_hapbams
        .map { meta, bam, bai, vcf, snvs, mods, hapbam, hapbai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            basecall_model: meta.basecall_model]
            return[new_meta, [[type: meta.type], hapbam], [[type: meta.type], hapbai]]
        }
        .groupTuple()
        .map{ meta, bam, bai ->
            def normal_bam = bam[0][0].type == "normal" ? bam[0][1] : bam[1][1]
            def tumor_bam = bam[0][0].type == "tumor" ? bam[0][1] : bam[1][1]
            def normal_bai = bai[0][0].type == "normal" ? bai[0][1] : bai[1][1]
            def tumor_bai = bai[0][0].type == "tumor" ? bai[0][1] : bai[1][1]
            // Return channel
            return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
        }
        .join(LONGPHASE_PHASE.out.vcf)
        .set{tumor_normal_severus}

    // Get ClairS input channel
    tumor_normal_severus
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
            def model = (!meta.clairS_model || meta.clairS_model.toString().trim() in ['', '[]']) ? clairs_modelMap.get(meta.basecall_model.toString().trim()) : meta.clairS_model
            return[meta , tumor_bam, tumor_bai, normal_bam, normal_bai, model]
        }
        .set { clairs_input }

    //
    // MODULE: CLAIRS
    //

    CLAIRS (
        clairs_input,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(CLAIRS.out.versions)
    // tumor_normal_severus -> meta:       [id, paired_data, platform, sex, fiber, basecall_model]
    //                         tumor_bam:  haplotagged aligned bam for tumor
    //                         tumor_bai:  indexes for tumor bam files
    //                         normal_bam: haplotagged aligned bam files for normal
    //                         normal_bai: indexes for normal bam files
    //                         phased_vcf: phased small variant vcf for normal

      // Get ClairS input channel

    tumor_normal_severus
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
            def model = (!meta.clairS_model || meta.clairS_model.toString().trim() in ['', '[]']) ? clairs_modelMap.get(meta.basecall_model.toString().trim()) : meta.clairS_model
            return[meta , tumor_bam, tumor_bai, normal_bam, normal_bai,model]
        }
        .set { clairs_input }

    //
    // MODULE: CLAIRS
    //

    CLAIRS (
        clairs_input,
        ch_fasta,
        ch_fai
    )

    SOMATIC_VEP (
        CLAIRS.out.vcf,
        params.genome,
        "homo_sapiens",
        111,
        '',
        fasta,
        []
    )

    ch_versions = ch_versions.mix(CLAIRS.out.versions)

    emit:
    tumor_normal_severus
    versions = ch_versions

}
