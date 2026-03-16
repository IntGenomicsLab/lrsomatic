include { CLAIRSTO                  } from '../../modules/local/clairsto/main.nf'
include { VCFSPLIT                  } from '../../modules/local/vcfsplit/main.nf'
include { LONGPHASE_PHASE           } from '../../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG        } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main.nf'
include { DEEPVARIANT               } from '../../subworkflows/nf-core/deepvariant/main.nf'
include { DEEPSOMATIC               } from '../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS   } from '../../subworkflows/local/small_variant_consensus.nf'


workflow TUMOR_ONLY_HAPPHASE {

    take:
    tumor_bams
    fasta
    fai
    dbsnp
    colors
    onekgenomes
    gnomad

    main:

    ch_versions = channel.empty()
    tumor_only_severus = channel.empty()
    somatic_vep = channel.empty()
    germline_vep = channel.empty()

    tumor_bams
        .map{ meta, bam, bai ->
            return [meta, bam, bai, meta.clairSTO_model]
        }
        .set{ tumor_bams }
    // [meta, bam, bai, clairSTO_model]  -- ClairS-TO model string appended for CLAIRSTO input

    tumor_bams
        .map { meta, bam, bai, _clairSTO_model ->
            def intervals = []
            return [meta,bam,bai, intervals]
        }
        .set{tumor_bams_deepvar}
    
    tumor_bams
        .map { meta, tumor_bam, tumor_bai, _clairSTO_model ->
            def normal_bam = []
            def normal_bai = []
            return [meta,normal_bam,normal_bai,tumor_bam,tumor_bai]
        }
        .set{tumor_bams_deepsomatic}

    DEEPVARIANT (
        tumor_bams_deepvar,
        fasta,
        fai,
        [[:],[]],
        [[:],[]]
    )

    DEEPSOMATIC (
        tumor_bams_deepsomatic,
        [[:],[]],
        fasta,
        fai,
        [[:],[]]
    )

    //
    // MODULE: CLAIRSTO
    //
    // call somatic/non-somatic variants
    // (* not called as germline * just non-somatic)

    CLAIRSTO (
        tumor_bams,
        fasta,
        fai,
        dbsnp,
        colors,
        onekgenomes,
        gnomad
    )

    CLAIRSTO.out.indel_vcf
                .join(CLAIRSTO.out.snv_vcf)
                .set{ clairsto_vcf }
    // [meta, indel_vcf, snv_vcf]  -- raw ClairS-TO variant calls

    //
    // MODULE: VCFSPLIT
    //
    // ClairSTO gives outputs in snv.vcf and indel.vcf
    // reformats them to be in somatic.vcf and nonsomatic.vcf

    VCFSPLIT (
        clairsto_vcf
    )

    VCFSPLIT.out.germline_vcf
        .join(VCFSPLIT.out.germline_tbi)
        .map { meta, vcf, tbi ->
            def new_meta = meta + [caller:'clairs-to']
            return [ new_meta, vcf, tbi]
        }
        .set{clairsto_germline_ch}
    
    DEEPVARIANT.out.vcf
        .join(DEEPVARIANT.out.vcf_index)
        .map{ meta, vcf, tbi ->
            def new_meta = meta + [caller:'deepvariant']
            return [new_meta, vcf, tbi]
        }
        .set{deepvariant_ch}

    clairsto_germline_ch
        .mix(deepvariant_ch)
        .set{mixed_vcfs}

    SMALL_VARIANT_CONSENSUS(
        mixed_vcfs,
        fasta,
        fai,
        params.germline_var_keep
    )
    // Add the nonsomatic vcf info
    // remove model info
    tumor_bams
        .join(SMALL_VARIANT_CONSENSUS.out.vcf)
        .map{ meta, bam, bai, _model, snps ->
            def svs = []
            def mods = []
            return[meta, bam, bai, snps, svs, mods]
        }
        .set{ tumor_bams_germlinevcf }
    // [meta, bam, bai, nonsomatic_vcf, [], []]  -- non-somatic variants used for phasing; svs and mods are empty placeholders for LONGPHASE_PHASE input

    VCFSPLIT.out.somatic_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta,vcf, extra]
        }
        .set { somatic_vep }
    // [meta, somatic_vcf, []]  -- PASS (somatic) variants for VEP annotation

    VCFSPLIT.out.germline_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta,vcf, extra]
        }
        .set { germline_vep }
    // [meta, germline_vcf, []]  -- non-somatic variants (relabelled PASS) for VEP annotation

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

    // Add phased nonsomatic vcf info
    // remove model info
    tumor_bams
        .join(LONGPHASE_PHASE.out.snv_vcf)
        .map { meta, bam, bai, _model, vcf ->
            def new_meta = meta + [type: "tumor"]
            def svs = []
            def mods = []
            return [new_meta, bam, bai, vcf, svs, mods]
        }
        .set{ tumor_bams_phasedvcf }
    // [meta+{type:"tumor"}, bam, bai, phased_nonsomatic_vcf, [], []]  -- type added; svs and mods are empty placeholders for LONGPHASE_HAPLOTAG

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

    // grab phased bams
    LONGPHASE_HAPLOTAG.out.bam
        .set{ haplotagged_bams }
    // [meta+{type:"tumor"}, haplotagged_bam]

    //
    // MODULES: SAMTOOLS_INDEX
    //
    // index the haplotagged bams
    SAMTOOLS_INDEX (
        haplotagged_bams
    )

    // join information and the phased VCF file
    haplotagged_bams
        .join(SAMTOOLS_INDEX.out.bai)
        .join(LONGPHASE_PHASE.out.snv_vcf)
        .join(LONGPHASE_PHASE.out.snv_vcf_index)
        .map{ meta, hap_bam, hap_bai, vcf, tbi ->
            def new_meta = [id: meta.id,
                            paired_data: meta.paired_data,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            clair3_model: meta.clair3_model,
                            clairS_model: meta.clairS_model,
                            clairSTO_model: meta.clairSTO_model,
                            kinetics: meta.kinetics]
            return [new_meta, hap_bam, hap_bai, [], [], vcf, tbi]
            }
        .set{ tumor_only_severus }
    // [meta, hap_bam, hap_bai, [], [], phased_vcf, phased_tbi]  -- normal_bam and normal_bai are [] (tumor-only mode)

    emit:
    tumor_only_severus
    somatic_vep
    germline_vep
    versions = ch_versions

}
