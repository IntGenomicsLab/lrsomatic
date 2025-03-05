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
    mixed_bams
    .branch{ meta, bam, bai ->
        normal: meta.type == "normal"
        tumor: meta.type == "tumor"
    }
    .set{mixed_bams}
    
    mixed_bams.normal
    .map{meta, bam, bai ->
        def new_meta = [id: meta.id,
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model]
        return[new_meta, bam, bai]
    }
    .set{normal_bams}

    mixed_bams.tumor
    .map{meta, bam, bai ->
        def new_meta = [id: meta.id,
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model]
        return[new_meta, bam, bai]
    }
    .set{tumor_bams}

    CLAIR3(
        normal_bams,
        fasta,
        fai
    )

    normal_bams
    .join(CLAIR3.out.germline_vcf)
    .map{ meta, bam, bai, vcf ->
        def svs = []
        def mods = []
        return [meta, bam, bai, vcf, svs, mods]
    }
    .set{normalbams_germlinevcf}

    LONGPHASE_PHASE(
        normalbams_germlinevcf,
        fasta,
        fai
    )
    
    normal_bams
    .join(LONGPHASE_PHASE.out.vcf)
    .map{meta, bam, bai, vcf ->
        def new_meta = meta + [type: "normal"]
        def snvs = []
        def mods = []
        return[new_meta, bam, bai, vcf, snvs, mods]
    }
    .set{normal_bams}

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

    LONGPHASE_HAPLOTAG(
        mixed_bams_vcf,
        fasta,
        fai
    )

    LONGPHASE_HAPLOTAG.out.bam
    .set{mixed_hapbams}

    SAMTOOLS_INDEX(
        mixed_hapbams
    )

    mixed_bams_vcf
    .join(mixed_hapbams)
    .join(SAMTOOLS_INDEX.out.bai)
    .set{mixed_hapbams}

    mixed_hapbams
    .map{ meta, bam, bai, vcf, snvs, mods, hapbam, hapbai->
        def new_meta = [id: meta.id,
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model]
        return[new_meta , [[type: meta.type], hapbam], [[type: meta.type], hapbai]]
    
    }
    .groupTuple()
    .map{meta, bam, bai ->

        def normal_bam = bam[0][0].type == "normal" ? bam[0][1] : bam[1][1]
        def tumor_bam = bam[0][0].type == "tumor" ? bam[0][1] : bam[1][1]
        def normal_bai = bai[0][0].type == "normal" ? bai[0][1] : bai[1][1]
        def tumor_bai = bai[0][0].type == "tumor" ? bai[0][1] : bai[1][1]

        return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    }
    .join(LONGPHASE_PHASE.out.vcf)
    .set{tumor_normal_severus}


    emit:
    tumor_normal_severus
    
}