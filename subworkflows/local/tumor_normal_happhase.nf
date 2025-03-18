include { CLAIR3 } from '../../modules/local/clair3/main.nf'
include { LONGPHASE_PHASE } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'

workflow TUMOR_NORMAL_HAPPHASE {
    take:
    mixed_bams
    fasta
    fai

    main:
    
    ch_versions = Channel.empty()
    
    // Branch input bams in normal and tumour
    mixed_bams
        .branch{ meta, bam, bai ->
            normal: meta.type == "normal"
            tumor: meta.type == "tumor"
        }
        .set{mixed_bams}
    
    // Get normal bams
    mixed_bams.normal
        .map{ meta, bam, bai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            basecall_model: meta.basecall_model]
            return[new_meta, bam, bai]
        }
        .set{normal_bams}
    
    // Get tumour bams
    mixed_bams.tumor
        .map{ meta, bam, bai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            basecall_model: meta.basecall_model]
            return[new_meta, bam, bai]
        }
        .set{tumor_bams}
    
    //
    // MODULE: CLAIR3
    //
    
    CLAIR3 (
        normal_bams,
        fasta,
        fai
    )
    
    ch_versions = ch_versions.mix(CLAIR3.out.versions)
    
    // Add germline vcf to normal bams
    normal_bams
        .join(CLAIR3.out.germline_vcf)
        .map{ meta, bam, bai, vcf ->
            def svs = []
            def mods = []
            return [meta, bam, bai, vcf, svs, mods]
        }
        .set{normalbams_germlinevcf}
    
    //
    // MODULE: LONGPHASE_PHASE
    //
    
    LONGPHASE_PHASE(
        normalbams_germlinevcf,
        fasta,
        fai
    )
    
    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)
    
    // Add phased vcf to normal bams
    normal_bams
        .join(LONGPHASE_PHASE.out.vcf)
        .map{meta, bam, bai, vcf ->
            def new_meta = meta + [type: "normal"]
            def snvs = []
            def mods = []
            return[new_meta, bam, bai, vcf, snvs, mods]
        }
        .set{normal_bams}
    
    // Add phased vcf to tumour bams
    tumor_bams
        .join(LONGPHASE_PHASE.out.vcf)
        .map{meta, bam, bai, vcf ->
            def new_meta = meta + [type: "tumor"]
            def snvs = []
            def mods = []
            return [new_meta, bam, bai, vcf, snvs, mods]
        }
        .mix(normal_bams)
        .set{mixed_bams_vcf}
    
    //
    // MODULE: LONGPHASE_HAPLOTAG
    //
    
    LONGPHASE_HAPLOTAG (
        mixed_bams_vcf,
        fasta,
        fai
    )
    
    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)
    
    // Get final tagged bams
    LONGPHASE_HAPLOTAG.out.bam
        .set{ mixed_hapbams }
    
    //
    // MODULE: SAMTOOLS_INDEX
    //
    
    SAMTOOLS_INDEX (
        mixed_hapbams
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    
    // Add index to channel
    mixed_bams_vcf
        .join(mixed_hapbams)
        .join(SAMTOOLS_INDEX.out.bai)
        .set{mixed_hapbams}
    
    // Group everything back together in one channel
    // Format of channel is [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf]
    mixed_hapbams
        .map { meta, bam, bai, vcf, snvs, mods, hapbam, hapbai ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            basecall_model: meta.basecall_model]
            return[new_meta , [[type: meta.type], hapbam], [[type: meta.type], hapbai]]
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

    emit:
    tumor_normal_severus
    versions = ch_versions
    
}