// Import modules
include { LONGPHASE_PHASE as LONGPHASE_PHASE_GERMLINE       } from '../../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG                                } from '../../../modules/nf-core/longphase/haplotag/main.nf'
include { LONGPHASE_MODCALL as LONGPHASE_MODCALL_GERMLINE   } from '../../../modules/local/longphase/modcall/main.nf'
include { SAMTOOLS_INDEX                                    } from '../../../modules/nf-core/samtools/index/main.nf'
include { BCFTOOLS_REHEADER as REHEADER_GERMLINE            } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_NORM as REJOIN_GERMLINE                  } from '../../../modules/nf-core/bcftools/norm/main'
include { PHASE_SOMATIC_FROM_READS                          } from '../../../modules/local/phase_somatic_from_reads/main.nf'


workflow PHASING_HAPLOTYPING {
    take:
    tumor_normal_bams // [meta, bam, bai]  -- all samples: tumor, normal, and tumor-only
    germline_vcf      // [meta, vcf, tbi]  -- germline small variants (from PAIRED_SMALLVAR_GERMLINE or TUMORONLY_SMALLVAR)
    somatic_vcf       // [meta, vcf, tbi]  -- somatic small variants (from PAIRED_SMALLVAR_SOMATIC or TUMORONLY_SMALLVAR)
    fasta             // [[:], fasta]
    fai               // [[:], fai]

    main:

    // Meta keys shared by every stream once 'type' is dropped, so BAMs, VCFs and modcall output join
    def meta_keys = ['id', 'paired_data', 'platform', 'sex', 'fiber', 'clair3_model', 'clairS_model', 'clairSTO_model']

    // SPLIT INTO PAIRED AND TUMOR ONLY
    // paired_data is set to the matched sample ID for paired samples, null/false for tumor-only
    tumor_normal_bams
        .branch { meta, _bams, _bai ->
            paired:      meta.paired_data
            tumor_only: !meta.paired_data
        }
        .set { branched_bams }
    // branched_bams.paired:     [meta, bam, bai]  -- tumor + normal from paired runs
    // branched_bams.tumor_only: [meta, bam, bai]  -- tumor-only samples

    // Strip 'type' from tumor-only meta (no type distinction needed in this stream)
    branched_bams.tumor_only
        .map { meta, bam, bai -> [ meta.subMap(meta_keys), bam, bai ] }
        .set{ tumor_only_ch }
    // tumor_only_ch: [meta (no type), bam, bai]

    // Split paired samples into normal and tumor streams for separate handling
    branched_bams.paired
        .branch { meta, _bam, _bai ->
            normal: meta.type == "normal"
            tumor:  meta.type == "tumor"
        }
        .set {paired_ch_branched}
    // paired_ch_branched.normal: [meta, bam, bai]  -- normal BAMs from T/N pairs
    // paired_ch_branched.tumor:  [meta, bam, bai]  -- tumor BAMs from T/N pairs

    // Strip 'type' from paired normal/tumor meta to allow joining with tumor-only channel
    paired_ch_branched.normal
        .map { meta, bam, bai -> [ meta.subMap(meta_keys), bam, bai ] }
        .set{ paired_normal_ch }
    paired_ch_branched.tumor
        .map { meta, bam, bai -> [ meta.subMap(meta_keys), bam, bai ] }
        .set{ paired_tumor_ch }
    // paired_normal_ch / paired_tumor_ch: [meta (no type), bam, bai]

    // Germline phasing uses normal BAMs (+ tumor-only BAMs used as their own "normal" proxy)
    tumor_only_ch
        .mix(paired_normal_ch)
        .set { normal_bams_w_tumoronly_ch }
    // normal_bams_w_tumoronly_ch: [meta, bam, bai]

    // MODCALL: Longphase base-modification calls, used as extra phasing evidence

    if (!params.skip_modcall) {

        //
        // MODULE: LONGPHASE_MODCALL_GERMLINE (label: process_high)
        // Input:  [meta, bam, bai]  -- normal BAMs (+ tumor-only BAMs)
        //         fasta / fai
        // Output: .mod_vcf -- [meta, vcf]  -- base modification calls (e.g. CpG methylation)
        //
        LONGPHASE_MODCALL_GERMLINE (
            normal_bams_w_tumoronly_ch,
            fasta,
            fai
        )
    }

    //
    // MODULE: BCFTOOLS_REHEADER (label: process_low), aliased REHEADER_GERMLINE
    // LongPhase haplotag only writes the contigs listed in the ##contig header of the phased germline VCF, so
    // reads on any other contig (chrM, unplaced/alt contigs, chrY without calls) are dropped from the haplotagged
    // BAM. ClairS-TO only lists contigs that have calls (chr1-22,X,Y by default), so replace the header contigs
    // with the full reference .fai. Records are untouched.
    // Input:  [meta, vcf, header, samples] / [meta2, fai]  (header and samples unused)
    // Output: .vcf + .index -- [meta, vcf.gz] / [meta, tbi]
    //
    REHEADER_GERMLINE (
        germline_vcf.map { meta, vcf, _tbi -> [ meta, vcf, [], [] ] },
        fai
    )
    REHEADER_GERMLINE.out.vcf
        .join(REHEADER_GERMLINE.out.index, failOnMismatch: true, failOnDuplicate: true)
        .set{ reheadered_germline_vcf }
    // reheadered_germline_vcf: [meta, vcf, tbi]  -- germline VCF with all reference contigs in the header

    //
    // MODULE: BCFTOOLS_NORM (label: process_medium), aliased REJOIN_GERMLINE
    // Join records sharing a POS (-m +any); LongPhase cannot take two records at one POS.
    //
    REJOIN_GERMLINE( reheadered_germline_vcf, fasta )

    REJOIN_GERMLINE.out.vcf
        .join(REJOIN_GERMLINE.out.index, failOnMismatch: true, failOnDuplicate: true)
        .set{ rejoined_germline_vcf }
    // rejoined_germline_vcf: [meta, vcf, tbi]

    // PHASING: germline (normal BAMs + germline VCF) builds the phase blocks every BAM is haplotagged with
    if (!params.skip_modcall) {
        normal_bams_w_tumoronly_ch
            .join(rejoined_germline_vcf)
            .join(LONGPHASE_MODCALL_GERMLINE.out.mod_vcf)
            .map { meta, bam, bai, vcf, _tbi, mods -> [ meta, bam, bai, vcf, [], mods ] }
            .set{ longphase_phase_germline_input_ch }
    }
    else {
        normal_bams_w_tumoronly_ch
            .join(rejoined_germline_vcf)
            .map { meta, bam, bai, vcf, _tbi -> [ meta, bam, bai, vcf, [], [] ] }
            .set{ longphase_phase_germline_input_ch }
    }
    // longphase_phase_germline_input_ch: [meta, bam, bai, germline_vcf, [] (svs), mod_vcf or []]

    //
    // MODULE: LONGPHASE_PHASE_GERMLINE (label: process_medium)
    // Input:  [meta, bam, bai, vcf, svs, mods]  -- normal BAMs + germline VCF (± mod VCF)
    //         fasta / fai
    // Output: .snv_vcf       -- [meta, vcf]  -- phased germline SNV VCF (PS tags added)
    //         .snv_vcf_index -- [meta, tbi]
    //
    LONGPHASE_PHASE_GERMLINE (
        longphase_phase_germline_input_ch,
        fasta,
        fai
    )

    LONGPHASE_PHASE_GERMLINE.out.snv_vcf
        .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf_index)
        .set{ phased_germline_vcf }
    // phased_germline_vcf: [meta, vcf, tbi]  -- Longphase-phased germline VCF

    // HAPLOTAGGING: every sample is tagged from the germline phase blocks; 'type' goes back into meta for output naming
    // A paired tumour and its normal share the type-less meta, so they meet the germline VCF via combine, not join
    tumor_only_ch.map { meta, bam, bai -> [ meta, bam, bai, 'tumor' ] }
        .mix(paired_tumor_ch.map { meta, bam, bai -> [ meta, bam, bai, 'tumor' ] })
        .mix(paired_normal_ch.map { meta, bam, bai -> [ meta, bam, bai, 'normal' ] })
        .combine(LONGPHASE_PHASE_GERMLINE.out.snv_vcf, by: 0)
        .set { haplotag_bams_vcf_ch }
    // haplotag_bams_vcf_ch: [meta (no type), bam, bai, type, phased_germline_vcf]

    if (!params.skip_modcall) {
        haplotag_bams_vcf_ch
            .combine(LONGPHASE_MODCALL_GERMLINE.out.mod_vcf.map { meta, mods -> [ meta.subMap(meta_keys), mods ] }, by: 0)
            .map { meta, bam, bai, type, vcf, mods -> [ meta + [ type: type ], bam, bai, vcf, [], mods ] }
            .set { longphase_haplotag_input_ch }
    }
    else {
        haplotag_bams_vcf_ch
            .map { meta, bam, bai, type, vcf -> [ meta + [ type: type ], bam, bai, vcf, [], [] ] }
            .set { longphase_haplotag_input_ch }
    }
    // longphase_haplotag_input_ch: [meta(+type), bam, bai, phased_germline_vcf, [], mod_vcf_or_[]]
    //   -- all samples (tumor-only, paired tumor, paired normal)

    //
    // MODULE: LONGPHASE_HAPLOTAG (label: process_medium)
    // Input:  [meta, bam, bai, phased_vcf, svs, mods]  -- BAM + phased germline VCF (± mod VCF)
    //         fasta / fai
    // Output: .bam -- [meta, bam]  -- BAM with HP (haplotype) and PS (phase set) tags added to reads
    //
    LONGPHASE_HAPLOTAG (
        longphase_haplotag_input_ch,
        fasta,
        fai
    )

    //
    // MODULE: SAMTOOLS_INDEX (label: process_medium)
    // Input:  [meta, bam]  -- haplotagged BAM
    // Output: .index -- [meta, index]  -- .bai or .crai, whichever matches the haplotag output format
    //
    SAMTOOLS_INDEX (
        LONGPHASE_HAPLOTAG.out.bam
    )
    LONGPHASE_HAPLOTAG.out.bam
        .join(SAMTOOLS_INDEX.out.index)
        .set{ tumor_normal_hapbams_ch }
    // tumor_normal_hapbams_ch: [meta, bam, bai]  -- haplotagged BAM with index

    //
    // MODULE: PHASE_SOMATIC_FROM_READS (label: process_low)
    // Somatic GT/PS from the HP/PS of each call's alt reads in the haplotagged tumour BAM, so they share the germline phase.
    // Input:  [meta, somatic_vcf, tbi, hap_bam, index] / fasta / fai
    // Output: .vcf / .tbi -- [meta, vcf.gz] / [meta, tbi]
    //
    // No failOnMismatch: a sample whose haplotagging failed (errorStrategy 'ignore') has no phased somatic VCF
    somatic_vcf
        .join(
            tumor_normal_hapbams_ch
                .filter { meta, _bam, _bai -> meta.type == 'tumor' }
                .map { meta, bam, bai -> [ meta.subMap(meta_keys), bam, bai ] },
            failOnDuplicate: true
        )
        .set{ somatic_phase_input_ch }
    // somatic_phase_input_ch: [meta, somatic_vcf, tbi, hap_bam, index]

    PHASE_SOMATIC_FROM_READS ( somatic_phase_input_ch, fasta, fai )

    PHASE_SOMATIC_FROM_READS.out.vcf
        .join(PHASE_SOMATIC_FROM_READS.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .set{ phased_somatic_vcf }
    // phased_somatic_vcf: [meta, vcf.gz, tbi]  -- somatic VCF phased from the haplotagged tumour reads


    emit:
    tumor_normal_hapbams_ch   // [meta, bam, bai]  -- haplotagged BAMs for all samples
    phased_germline_vcf       // [meta, vcf, tbi]  -- phased germline VCF (used by SEVERUS + VEP)
    phased_somatic_vcf        // [meta, vcf, tbi]  -- phased somatic VCF (used by VEP)
}
