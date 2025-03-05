include { CLAIRSTO } from '../../modules/local/clairsto/main.nf'
include { VCFSPLIT } from '../../modules/local/vcfsplit/main.nf'
include { LONGPHASE_PHASE          } from '../../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'

workflow TUMOR_ONLY_HAPPHASE {

    take:
    tumor_bams
    fasta
    fai

    main:
    tumor_bams.view()
    CLAIRSTO(
        tumor_bams,
        fasta,
        fai
    )

    VCFSPLIT(
        CLAIRSTO.out.vcfs
    )

    tumor_bams
    .join(VCFSPLIT.out.germline_vcf)
    .map{ meta, bam, bai, snps ->
        def svs = []
        def mods = []
        return[meta, bam, bai, snps, svs, mods]
    }
    .set{tumor_bams_germlinevcf}
    // [meta, bam, bai, vcf]

    LONGPHASE_PHASE(
        tumor_bams_germlinevcf,
        fasta,
        fai
    )

    tumor_bams
    .join(LONGPHASE_PHASE.out.vcf)
    .map{ meta, bam, bai, snps ->
        def svs = []
        def mods = []
        return [meta, bam, bai, snps, svs, mods]
    }
    .set{tumor_bams_phasedvcf}

    LONGPHASE_HAPLOTAG(
        tumor_bams_phasedvcf,
        fasta,
        fai
    )

    LONGPHASE_HAPLOTAG.out.bam
        .set{haplotagged_bams}

    SAMTOOLS_INDEX(
        haplotagged_bams
    )

    haplotagged_bams
    .join(SAMTOOLS_INDEX.out.bai)
    .join(LONGPHASE_PHASE.out.vcf)
    .map{meta, hap_bam, hap_bai, vcf ->
        def new_meta = [id: meta.id,
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model]
        return[new_meta, hap_bam, hap_bai, [],[], vcf]
        }
    .set{tumor_only_severus}

    emit:
    tumor_only_severus
    
}