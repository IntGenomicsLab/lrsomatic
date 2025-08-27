include { CLAIRSTO                  } from '../../modules/local/clairsto/main.nf'
include { VCFSPLIT                  } from '../../modules/local/vcfsplit/main.nf'
include { LONGPHASE_PHASE           } from '../../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG        } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'
include {ENSEMBLVEP_VEP as SOMATIC_VEP} from '../../modules/nf-core/ensemblvep/vep/main.nf'
include {ENSEMBLVEP_VEP as GERMLINE_VEP} from '../../modules/nf-core/ensemblvep/vep/main.nf'

workflow TUMOR_ONLY_HAPPHASE {

    take:
    tumor_bams
    fasta
    fai
    clairSTO_modelMap

    main:

    ch_versions = Channel.empty()

    ch_versions = Channel.empty()

    tumor_bams
        .map{ meta, bam, bai ->
           def clairSTO_model = (!meta.clairSTO_model || meta.clairSTO_model.toString().trim() in ['', '[]']) ? clairSTO_modelMap.get(meta.basecall_model.toString().trim()) : meta.clairSTO_model
            return [meta, bam, bai, clairSTO_model]
        }
        .set{ tumor_bams }

    //
    // MODULE: CLAIRSTO
    //
    // call somatic/non-somatic variants
    // (* not called as germline * just non-somatic)

    CLAIRSTO (
        tumor_bams,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(CLAIRSTO.out.versions)

    CLAIRSTO.out.indel_vcf
                .join(CLAIRSTO.out.snv_vcf)
                .set{ clairsto_vcf }

    ch_versions = ch_versions.mix(CLAIRSTO.out.versions)
    // clairsto_vcf -> meta:      [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                 indel_vcf: vcf for indels
    //                 snv_vcf:   vcf for snvs

    //
    // MODULE: VCFSPLIT
    //
    // ClairSTO gives outputs in snv.vcf and indel.vcf
    // reformats them to be in somatic.vcf and nonsomatic.vcf

    VCFSPLIT (
        clairsto_vcf
    )
    ch_versions = ch_versions.mix(VCFSPLIT.out.versions)

    ch_versions = ch_versions.mix(VCFSPLIT.out.versions)

    // Add the nonsomatic vcf info
    // remove model info
    tumor_bams
        .join(VCFSPLIT.out.germline_vcf)
        .map{ meta, bam, bai, model, snps ->
            def svs = []
            def mods = []
            return[meta, bam, bai, snps, svs, mods]
        }
        .set{ tumor_bams_germlinevcf }


    VCFSPLIT.out.somatic_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta,vcf, extra]
        }
        .set { somatic_vep }

    VCFSPLIT.out.germline_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta,vcf, extra]
        }
        .set { germline_vep }

    SOMATIC_VEP (
        somatic_vep,
        params.genome,
        "homo_sapiens",
        111,
        '',
        fasta,
        []
    )

    GERMLINE_VEP (
        germline_vep,
        params.genome,
        "homo_sapiens",
        111,
        '',
        fasta,
        []
    )


    // tumor_bams_germlinevcf -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                           bam:  list of concatenated aligned bams
    //                           bai:  indexes for bam files
    //                           vcf:  tumor small nonsomatic variant vcf
    //                           svs:  structural variant vcf (empty)
    //                           mods: modcall-generated VCF with modifications (empty)

    //
    // MODULES: LONGPHASE_PHASE
    //
    // Phase tumor bams on nonsomatic vcf

    LONGPHASE_PHASE (
        tumor_bams_germlinevcf,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

    // Add phased nonsomatic vcf info
    // remove model info
    tumor_bams
        .join(LONGPHASE_PHASE.out.vcf)
        .map{ meta, bam, bai, model, snps ->
            def svs = []
            def mods = []
            return [meta, bam, bai, snps, svs, mods]
        }
        .set{ tumor_bams_phasedvcf }

    // tumor_bams_germlinevcf -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                           bam:  list of concatenated aligned bams
    //                           bai:  indexes for bam files
    //                           vcf:  phased tumor small nonsomatic variant vcf
    //                           svs:  structural variant vcf (empty)
    //                           mods: modcall-generated VCF with modifications (empty)

    //
    // MODULES: LONGPHASE_HAPLOTAG
    //
    // Haplotag the tumor bams

    LONGPHASE_HAPLOTAG (
        tumor_bams_phasedvcf,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    // grab phased bams
    LONGPHASE_HAPLOTAG.out.bam
        .set{ haplotagged_bams }

    // haplotagged_bams -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                     bams: list of concatenated aligned bams

    //
    // MODULES: SAMTOOLS_INDEX
    //
    // index the haplotagged bams
    SAMTOOLS_INDEX (
        haplotagged_bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // join information and the phased VCF file
    haplotagged_bams
        .join(SAMTOOLS_INDEX.out.bai)
        .join(LONGPHASE_PHASE.out.vcf)
        .map{meta, hap_bam, hap_bai, vcf ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            basecall_model: meta.basecall_model]
            return[new_meta, hap_bam, hap_bai, [],[], vcf]
            }
        .set{ tumor_only_severus }

    // tumor_normal_severus -> meta:       [id, paired_data, platform, sex, fiber, basecall_model]
    //                         hap_bam:  haplotagged aligned bam for tumor
    //                         hap_bai:  indexes for tumor bam files
    //                         normal_bam: haplotagged aligned bam files for normal (empty)
    //                         normal_bai: indexes for normal bam files (empty)
    //                         phased_vcf: phased small variant vcf

    emit:
    tumor_only_severus
    versions = ch_versions

}
