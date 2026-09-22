// Import modules
include { LONGPHASE_PHASE as LONGPHASE_PHASE_GERMLINE       } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_PHASE as LONGPHASE_PHASE_SOMATIC        } from '../../modules/nf-core/longphase/phase/main.nf'
include { LONGPHASE_HAPLOTAG                                } from '../../modules/nf-core/longphase/haplotag/main.nf'
include { LONGPHASE_MODCALL as LONGPHASE_MODCALL_GERMLINE   } from '../../modules/local/longphase/modcall/main.nf'
include { LONGPHASE_MODCALL as LONGPHASE_MODCALL_SOMATIC    } from '../../modules/local/longphase/modcall/main.nf'
include { SAMTOOLS_INDEX                                    } from '../../modules/nf-core/samtools/index/main.nf'
include { BCFTOOLS_CONCAT                                   } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                     } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW                                     } from '../../modules/local/bcftools/view/main.nf'
include { VCFTAG as TAG_SOMATIC                             } from '../../modules/local/vcftag/main.nf'
include { VCFTAG as TAG_GERMLINE                            } from '../../modules/local/vcftag/main.nf'


workflow PHASING_HAPLOTYPING {
    take:
    tumor_normal_bams // [meta, bam, bai]  -- all samples: tumor, normal, and tumor-only
    germline_vcf      // [meta, vcf, tbi]  -- germline small variants (from PAIRED_SMALLVAR_GERMLINE or TUMORONLY_SMALLVAR)
    somatic_vcf       // [meta, vcf, tbi]  -- somatic small variants (from PAIRED_SMALLVAR_SOMATIC or TUMORONLY_SMALLVAR)
    fasta             // [[:], fasta]
    fai               // [[:], fai]

    main:

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

    branched_bams.paired
        .set{ paired_ch }

    // Strip 'type' from tumor-only meta (no type distinction needed in this stream)
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
    // tumor_only_ch: [meta (no type), bam, bai]

    // Split paired samples into normal and tumor streams for separate handling
    paired_ch
        .branch { meta, _bam, _bai ->
            normal: meta.type == "normal"
            tumor:  meta.type == "tumor"
        }
        .set {paired_ch_branched}
    // paired_ch_branched.normal: [meta, bam, bai]  -- normal BAMs from T/N pairs
    // paired_ch_branched.tumor:  [meta, bam, bai]  -- tumor BAMs from T/N pairs

    // Strip 'type' from paired normal/tumor meta to allow joining with tumor-only channel
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
    // paired_normal_ch: [meta (no type), bam, bai]

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
    // paired_tumor_ch: [meta (no type), bam, bai]

    // Germline phasing uses normal BAMs (+ tumor-only BAMs used as their own "normal" proxy)
    tumor_only_ch
        .mix(paired_normal_ch)
        .set { normal_bams_w_tumoronly_ch }
    // normal_bams_w_tumoronly_ch: [meta, bam, bai]
    //   -- normal BAMs from T/N pairs + tumor-only BAMs (both phased with germline VCF)

    // Somatic phasing uses tumor BAMs (+ tumor-only BAMs)
    tumor_only_ch
        .mix(paired_tumor_ch)
        .set{ tumor_bams_ch}
    // tumor_bams_ch: [meta, bam, bai]  -- tumor BAMs from T/N pairs + tumor-only BAMs

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

        //
        // MODULE: LONGPHASE_MODCALL_SOMATIC (label: process_high)
        // Input:  [meta, bam, bai]  -- tumor BAMs (+ tumor-only BAMs)
        //         fasta / fai
        // Output: .mod_vcf -- [meta, vcf]  -- base modification calls for tumor
        //

        LONGPHASE_MODCALL_SOMATIC (
            tumor_bams_ch,
            fasta,
            fai
        )

    }

    //
    // MODULE: VCFTAG (label: process_single), aliased TAG_SOMATIC / TAG_GERMLINE
    // Stamp each arm with an INFO provenance flag before the merge. This is the only point where
    // germline-vs-somatic origin is unambiguous for every caller: GERMLINE_CONSENSUS can emit
    // records that never passed through VCFSPLIT, so tagging earlier would leave holes. After the
    // merge the two populations are otherwise indistinguishable -- both carry FILTER=PASS.
    // LongPhase preserves custom INFO keys, so the flags survive phasing (verified on v2.0.1).
    //
    TAG_SOMATIC ( somatic_vcf,  'SOMATIC'  )
    TAG_GERMLINE( germline_vcf, 'GERMLINE' )

    TAG_SOMATIC.out.vcf
        .join(TAG_SOMATIC.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .set{ tagged_somatic_vcf }
    TAG_GERMLINE.out.vcf
        .join(TAG_GERMLINE.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .set{ tagged_germline_vcf }
    // tagged_*_vcf: [meta, vcf, tbi]

    // Somatic phasing needs germline and somatic sites in one VCF for consistent phase blocks
    tagged_germline_vcf
        .join(tagged_somatic_vcf)
        .map { meta, germ_vcf, germ_tbi, som_vcf, som_tbi ->
                def vcfs = [som_vcf, germ_vcf]  // somatic first (higher priority in phasing)
                def tbis = [som_tbi, germ_tbi]
                return [ meta, vcfs, tbis]
        }
        .set{germline_somatic_vcfs}
    // germline_somatic_vcfs (pre-concat): [meta, [somatic_vcf, germline_vcf], [somatic_tbi, germline_tbi]]

    //
    // MODULE: BCFTOOLS_CONCAT (label: process_medium)
    // Input:  [meta, [vcfs...], [tbis...]]  -- somatic + germline VCFs to concatenate
    // Output: .vcf -- [meta, vcf]  -- unsorted concatenated VCF
    //
        BCFTOOLS_CONCAT(germline_somatic_vcfs)
        BCFTOOLS_CONCAT.out.vcf
            .set{concat_out}
    // concat_out: [meta, vcf]  -- concatenated (unsorted) somatic+germline VCF

    //
    // MODULE: BCFTOOLS_SORT (label: process_medium)
    // Input:  [meta, vcf]  -- unsorted concatenated VCF
    // Output: .vcf -- [meta, vcf]  -- coordinate-sorted VCF
    //         .tbi -- [meta, tbi]
    //
        BCFTOOLS_SORT(concat_out)
        BCFTOOLS_SORT.out.vcf
            .set{germline_somatic_vcfs}
    // germline_somatic_vcfs (final): [meta, vcf]  -- sorted combined somatic+germline VCF for somatic phasing

    // PHASING: germline (normal BAMs + germline VCF) builds the phase blocks; somatic (tumor BAMs + merged VCF) transfers them
    if (!params.skip_modcall) {
        // With modcall: include base-modification VCF as additional phasing evidence
        normal_bams_w_tumoronly_ch
            .join(tagged_germline_vcf)
            .join(LONGPHASE_MODCALL_GERMLINE.out.mod_vcf)
            .map { meta, bam, bai, vcf, _tbi, mods->
                def svs = []  // SVs for phasing are not used here
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_germline_input_ch }
        // longphase_phase_germline_input_ch: [meta, bam, bai, germline_vcf, [], mod_vcf]

        tumor_bams_ch
            .join(germline_somatic_vcfs)
            .join(LONGPHASE_MODCALL_SOMATIC.out.mod_vcf)
            .map { meta, bam, bai, vcf, mods->
                def svs = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_somatic_input_ch }
        // longphase_phase_somatic_input_ch: [meta, bam, bai, somatic+germline_vcf, [], mod_vcf]
    }
    else {
        // Without modcall: empty lists for SVs and mods
        normal_bams_w_tumoronly_ch
            .join(tagged_germline_vcf)
            .map { meta, bam, bai, vcf, _tbi ->
                def svs = []
                def mods = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_germline_input_ch }
        // longphase_phase_germline_input_ch: [meta, bam, bai, germline_vcf, [], []]

        tumor_bams_ch
            .join(germline_somatic_vcfs)
            .map { meta, bam, bai, vcf ->
                def svs = []
                def mods = []
                return [ meta, bam, bai, vcf, svs, mods ]
            }
            .set{ longphase_phase_somatic_input_ch }
        // longphase_phase_somatic_input_ch: [meta, bam, bai, somatic+germline_vcf, [], []]
    }

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

    //
    // MODULE: LONGPHASE_PHASE_SOMATIC (label: process_medium)
    // Input:  [meta, bam, bai, combined_vcf, svs, mods]  -- tumor BAMs + somatic+germline VCF (± mod VCF)
    //         fasta / fai
    // Output: .snv_vcf       -- [meta, vcf]  -- phased somatic (+ germline) VCF
    //         .snv_vcf_index -- [meta, tbi]
    //
    LONGPHASE_PHASE_SOMATIC (
        longphase_phase_somatic_input_ch,
        fasta,
        fai
    )

    LONGPHASE_PHASE_SOMATIC.out.snv_vcf
        .join(LONGPHASE_PHASE_SOMATIC.out.snv_vcf_index)
        .set{ phased_somatic_germline_vcf }
    // phased_somatic_germline_vcf: [meta, vcf, tbi]  -- Longphase-phased somatic+germline VCF (unfiltered)

    //
    // MODULE: BCFTOOLS_VIEW (label: process_medium)
    // Reduce the phased somatic+germline VCF to the somatic arm, selecting on the INFO/SOMATIC flag
    // stamped before the merge -- by provenance, not position. The previous `-T <somatic vcf>`
    // targets file matched CHROM/POS only, so germline records co-located with a somatic call were
    // retained and became indistinguishable downstream. PS/HP tags on somatic variants survive;
    // germline records are dropped here but stay published under variants/phased/ and vep/germline/.
    // Input:  [meta, phased_combined_vcf, phased_combined_tbi]
    // Output: .vcf -- [meta, vcf.gz]  -- phased somatic-only VCF
    //         .tbi -- [meta, tbi]
    //
    BCFTOOLS_VIEW ( phased_somatic_germline_vcf )

    BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)
        .set{ phased_somatic_vcf }
    // phased_somatic_vcf: [meta, vcf.gz, tbi]  -- phased somatic-only VCF (germline removed)

    // HAPLOTAGGING: every sample is tagged from the germline phase blocks; 'type' goes back into meta for output naming

    if(!params.skip_modcall) {
        // Strip 'type' so modcall output joins the type-less channels
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
        // modcall_vcf_ch: [meta (no type), mod_vcf]  -- base modification VCF from germline modcall

            // Build haplotag input for tumor-only samples (re-add type:"tumor")
            tumor_only_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ tumor_only_ch }
            // tumor_only_ch (updated): [meta+type:tumor, bam, bai, phased_germline_vcf, [], mod_vcf]

            paired_tumor_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_tumor_ch }
            // paired_tumor_ch (updated): [meta+type:tumor, bam, bai, phased_germline_vcf, [], mod_vcf]

            paired_normal_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .join(modcall_vcf_ch)
                .map { meta, bam, bai, vcf, mods ->
                    def new_meta = meta + [type : "normal"]
                    def svs = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_normal_ch }
            // paired_normal_ch (updated): [meta+type:normal, bam, bai, phased_germline_vcf, [], mod_vcf]

    }
    else {
        // Without modcall: empty lists for mods
            tumor_only_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ tumor_only_ch }
            // tumor_only_ch (updated): [meta+type:tumor, bam, bai, phased_germline_vcf, [], []]

            paired_tumor_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "tumor"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_tumor_ch }
            // paired_tumor_ch (updated): [meta+type:tumor, bam, bai, phased_germline_vcf, [], []]

            paired_normal_ch
                .join(LONGPHASE_PHASE_GERMLINE.out.snv_vcf)
                .map { meta, bam, bai, vcf ->
                    def new_meta = meta + [type : "normal"]
                    def svs = []
                    def mods = []
                    return [new_meta, bam, bai, vcf, svs, mods]
                }
                .set{ paired_normal_ch }
            // paired_normal_ch (updated): [meta+type:normal, bam, bai, phased_germline_vcf, [], []]

    }

    // Merge all sample types for haplotagging in a single LONGPHASE_HAPLOTAG call
    tumor_only_ch
        .mix(paired_tumor_ch)
        .mix(paired_normal_ch)
        .set {longphase_haplotag_input_ch}
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

    LONGPHASE_HAPLOTAG.out.bam
        .set{ tumor_normal_hapbams_ch }
    // tumor_normal_hapbams_ch (pre-index): [meta, bam]  -- haplotagged BAM (no index yet)

    //
    // MODULE: SAMTOOLS_INDEX (label: process_medium)
    // Input:  [meta, bam]  -- haplotagged BAM
    // Output: .bai -- [meta, bai]
    //
    SAMTOOLS_INDEX (
        tumor_normal_hapbams_ch
    )
    tumor_normal_hapbams_ch
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ tumor_normal_hapbams_ch }
    // tumor_normal_hapbams_ch (final): [meta, bam, bai]  -- haplotagged BAM with index


    emit:
    tumor_normal_hapbams_ch   // [meta, bam, bai]  -- haplotagged BAMs for all samples
    phased_germline_vcf       // [meta, vcf, tbi]  -- phased germline VCF (used by SEVERUS + VEP)
    phased_somatic_vcf        // [meta, vcf, tbi]  -- phased somatic VCF (used by VEP)
}
