// IMPORT MODULES
include { CLAIRSTO                  } from '../../../modules/local/clairsto/main.nf'
include { VCFSPLIT                  } from '../../../modules/local/vcfsplit/main.nf'

// IMPORT SUBWORKFLOWS
include { DEEPVARIANT                                   } from '../../../subworkflows/nf-core/deepvariant/main.nf'
include { DEEPSOMATIC                                   } from '../../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS as GERMLINE_CONSENSUS } from '../../../subworkflows/local/small_variant_consensus.nf'
include { SMALL_VARIANT_CONSENSUS as SOMATIC_CONSENSUS  } from '../../../subworkflows/local/small_variant_consensus.nf'


workflow TUMORONLY_SMALLVAR {

    take:
    tumor_bams           // [meta, tumor_bam, tumor_bai]  -- tumor-only aligned BAMs (no matched normal)
    fasta                // [[:], fasta]
    fai                  // [[:], fai]
    clairsto_pon_channel // [ [pon_vcf_path, ...], [is_population_allele_flag, ...] ]
    //                       used by ClairS-TO to filter germline variants with population allele databases
    ds_pon_channel       // [ [pon_vcf_path, ...] ] or [ [] ]
    //                       user-supplied DeepSomatic PON VCFs; empty list => container defaults

    main:

    somatic_vcf = channel.empty()
    germline_vcf = channel.empty()

    // CLAIRS-TO: somatic AND germline variant calling from tumor-only BAM
    // ClairS-TO uses a panel-of-normals / population allele database to separate somatic from germline
    // Runs if either somatic or germline clair calling is requested (produces both jointly)

    if(params.somatic_var_keep.contains('clair') || params.germline_var_keep.contains('clair')) {
        // Append model name and PoN info to build the full CLAIRSTO input
        tumor_bams
            .map { meta, bam, bai ->
                return [ meta, bam, bai, meta.clairSTO_model]
            }
            .combine(clairsto_pon_channel)
            .set{ clairsto_input_ch}
        // clairsto_input_ch: [meta, bam, bai, clairSTO_model_str, [pon_vcf_paths], [pon_flags]]

        //
        // MODULE: CLAIRSTO (label: process_high)
        // Input:  [meta, bam, bai, model_str, [pon_vcfs], [pon_flags]]
        //         fasta / fai
        // Output: .snv_vcf   -- [meta, vcf]  -- SNV calls (germline + somatic, unsplit)
        //         .indel_vcf -- [meta, vcf]  -- indel calls (germline + somatic, unsplit)
        //
        CLAIRSTO (
            clairsto_input_ch,
            fasta,
            fai
        )

        // SPLIT CLAIRSTO GERMLINE AND SOMATIC VARIATION
        // ClairS-TO outputs a combined VCF with FILTER tags indicating somatic/germline status;
        // VCFSPLIT separates these into two VCFs

        CLAIRSTO.out.indel_vcf
                    .join(CLAIRSTO.out.snv_vcf)
                    .set{ clairsto_combined_vcf }
        // clairsto_combined_vcf: [meta, indel_vcf, snv_vcf]

        //
        // MODULE: VCFSPLIT (label: process_single)
        // Input:  [meta, indel_vcf, snv_vcf]  -- combined ClairS-TO output
        // Output: .germline_vcf -- [meta, vcf]  -- germline variants only
        //         .germline_tbi -- [meta, tbi]
        //         .somatic_vcf  -- [meta, vcf]  -- somatic variants only
        //         .somatic_tbi  -- [meta, tbi]
        //
        VCFSPLIT (
            clairsto_combined_vcf
        )

        VCFSPLIT.out.germline_vcf
            .join(VCFSPLIT.out.germline_tbi)
            .map { meta, vcf, tbi ->
                def new_meta = meta + [caller:'clairs-to']
                return [ new_meta, vcf, tbi]
            }
            .set{clairsto_germline_ch}
        // clairsto_germline_ch: [meta(+caller:'clairs-to'), vcf, tbi]  -- germline variants

        VCFSPLIT.out.somatic_vcf
            .join(VCFSPLIT.out.somatic_tbi)
            .map { meta, vcf, tbi ->
                def new_meta = meta + [caller:'clairs-to']
                return [ new_meta, vcf, tbi]
            }
            .set{clairsto_somatic_ch}
        // clairsto_somatic_ch: [meta(+caller:'clairs-to'), vcf, tbi]  -- somatic variants
    }

    // DEEPVARIANT: germline-only variant calling (no somatic mode for tumor-only)
    if(params.germline_var_keep.contains('deepvariant')) {

        //
        // SUBWORKFLOW: DEEPVARIANT (nf-core)
        // Input:  [meta, bam, bai, []]  -- [] = genome-wide (no interval list)
        //         fasta / fai / [[:],[]] x2  -- empty PAR/GFF
        // Output: .vcf       -- [meta, vcf]
        //         .vcf_index -- [meta, tbi]
        //
        tumor_bams
            .map { meta, bam, bai  ->
                def intervals = []
                return [meta,bam,bai, intervals]
            }
            .set{deepvariant_input_ch}
        // deepvariant_input_ch: [meta, bam, bai, []]

        DEEPVARIANT (
            deepvariant_input_ch,
            fasta,
            fai,
            [[:],[]],  // PAR regions (not used)
            [[:],[]]   // GFF annotation (not used)
        )

        DEEPVARIANT.out.vcf
            .join(DEEPVARIANT.out.vcf_index)
            .map{ meta, vcf, tbi ->
                def new_meta = meta + [caller:'deepvariant']
                return [new_meta, vcf, tbi]
            }
            .set{deepvariant_ch}
        // deepvariant_ch: [meta(+caller:'deepvariant'), vcf, tbi]
    }

    // COMBINE GERMLINE VARIANTS
    // If both callers requested: run consensus; otherwise pass through single-caller output
    if (params.germline_var_keep.size() > 1) {
        clairsto_germline_ch
            .mix(deepvariant_ch)
            .set{combined_germline_ch}
        // combined_germline_ch: [meta(+caller), vcf, tbi]  -- one item per caller per sample

        // SUBWORKFLOW: GERMLINE_CONSENSUS (SMALL_VARIANT_CONSENSUS alias)
        GERMLINE_CONSENSUS(
            combined_germline_ch,
            fasta,
            fai,
            params.prioritize_caller_germline,
            params.germline_var_combine
        )
        GERMLINE_CONSENSUS.out.vcf
            .join(GERMLINE_CONSENSUS.out.tbi)
            .set{germline_vcf}
        // germline_vcf: [meta(+caller from consensus), vcf, tbi]
    }
    else if (params.germline_var_keep == ['clair']) {
        clairsto_germline_ch
            .set{germline_vcf}
    }
    else if (params.germline_var_keep == ['deepvariant']) {
        deepvariant_ch
            .set{germline_vcf}
    }

    // DEEPSOMATIC: somatic variant calling in tumor-only mode (no matched normal)
    // Normal BAM/BAI are passed as empty lists; DeepSomatic uses the model's internal normal baseline
    if(params.somatic_var_keep.contains('deepsomatic')) {
        tumor_bams
            .map { meta, tumor_bam, tumor_bai ->
                def normal_bam = []
                def normal_bai = []
                return [meta,normal_bam,normal_bai,tumor_bam,tumor_bai]
            }
            .set{deepsomatic_input_ch}
        // deepsomatic_input_ch: [meta, [], [], tumor_bam, tumor_bai]
        //   empty normal_bam/bai signals tumor-only mode to DEEPSOMATIC subworkflow

        //
        // SUBWORKFLOW: DEEPSOMATIC (local)
        // Input:  [meta, [], [], tumor_bam, tumor_bai]  -- tumor-only (no normal)
        //         [[:],[]] / fasta / fai / [[:],[]]
        // Output: .vcf       -- [meta, vcf]
        //         .vcf_index -- [meta, tbi]
        //
        DEEPSOMATIC (
            deepsomatic_input_ch,
            [[:],[]],  // intervals (empty = genome-wide)
            fasta,
            fai,
            [[:],[]],  // GZI (empty if FASTA is uncompressed)
            ds_pon_channel
        )
        DEEPSOMATIC.out.vcf
            .join(DEEPSOMATIC.out.vcf_index)
            .map{ meta, vcf, tbi ->
                def new_meta = meta + [caller:'deepsomatic']
                return [new_meta, vcf, tbi]
            }
            .set{deepsomatic_ch}
        // deepsomatic_ch: [meta(+caller:'deepsomatic'), vcf, tbi]
    }

    // COMBINE SOMATIC VARIATION
    if (params.somatic_var_keep.size() > 1) {
        clairsto_somatic_ch
            .mix(deepsomatic_ch)
            .set{combined_somatic_ch}
        // combined_somatic_ch: [meta(+caller), vcf, tbi]  -- one item per caller per sample

        // SUBWORKFLOW: SOMATIC_CONSENSUS (SMALL_VARIANT_CONSENSUS alias)
        SOMATIC_CONSENSUS(
            combined_somatic_ch,
            fasta,
            fai,
            params.prioritize_caller_somatic,
            params.somatic_var_combine
        )
        SOMATIC_CONSENSUS.out.vcf
            .join(SOMATIC_CONSENSUS.out.tbi)
            .set{somatic_vcf}
        // somatic_vcf: [meta(+caller from consensus), vcf, tbi]
    }
    else if (params.somatic_var_keep == ['clair']) {
        clairsto_somatic_ch
            .set{somatic_vcf}
    }
    else if (params.somatic_var_keep == ['deepsomatic']) {
        deepsomatic_ch
            .set{somatic_vcf}
    }

    // Strip 'caller' from meta before emitting both VCFs
    somatic_vcf
        .map{ meta, vcf, tbi  ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, vcf, tbi]
        }
        .set{somatic_vcf}

    germline_vcf
        .map{ meta, vcf, tbi  ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, vcf, tbi]
        }
        .set{germline_vcf}

    emit:
    somatic_vcf  // [meta, vcf, tbi]  -- final somatic VCF (ClairS-TO, DeepSomatic, or consensus)
    germline_vcf // [meta, vcf, tbi]  -- final germline VCF (ClairS-TO germline, DeepVariant, or consensus)


}
