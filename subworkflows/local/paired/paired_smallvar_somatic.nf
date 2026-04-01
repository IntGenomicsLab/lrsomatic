// IMPORT MODULES
include { CLAIRS                    } from '../../../modules/local/clairs/main.nf'
include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort'

// IMPORT SUBWORKFLOWS
include { DEEPSOMATIC                                    } from '../../../subworkflows/local/deepsomatic.nf'
include { SMALL_VARIANT_CONSENSUS as SOMATIC_CONSENSUS   } from '../../../subworkflows/local/small_variant_consensus.nf'

workflow PAIRED_SMALLVAR_SOMATIC {

    take:
    tumor_normal_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    fasta             // [[:], fasta]
    fai               // [[:], fai]

    main:
    somatic_vcf = channel.empty()

    // CLAIRS: somatic SNV/indel calling from T/N paired BAMs
    if(params.somatic_var_keep.contains('clair')) {
        // Append ClairS model name (from meta) as the last element for CLAIRS module
        tumor_normal_bams
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                return[meta , tumor_bam, tumor_bai, normal_bam, normal_bai, meta.clairS_model]
            }
            .set { clairs_input }
        // clairs_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, clairS_model_str]

        //
        // MODULE: CLAIRS (label: process_high)
        // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, model_str]
        //         fasta / fai
        // Output: .vcfs -- [meta, [snv_vcf, indel_vcf]]  -- separate SNV and indel VCFs
        //         .tbi  -- [meta, [snv_tbi, indel_tbi]]
        //
        CLAIRS (
            clairs_input,
            fasta,
            fai
        )

        // CONCAT CLAIRS INDEL AND SNV OUTPUT
        // ClairS outputs separate SNV and indel VCFs; merge into a single sorted VCF
        CLAIRS.out.vcfs
            .join(CLAIRS.out.tbi)
            .set{clairs_out}
        // clairs_out: [meta, [snv_vcf, indel_vcf], [snv_tbi, indel_tbi]]

        //
        // MODULE: BCFTOOLS_CONCAT (label: process_medium)
        // Input:  [meta, [vcf...], [tbi...]]
        // Output: .vcf -- [meta, vcf]  -- unsorted concatenated SNV+indel VCF
        //
        BCFTOOLS_CONCAT (
            clairs_out
        )

        //
        // MODULE: BCFTOOLS_SORT (label: process_medium)
        // Input:  [meta, vcf]
        // Output: .vcf -- [meta, vcf]  -- coordinate-sorted VCF
        //         .tbi -- [meta, tbi]
        //
        BCFTOOLS_SORT (
            BCFTOOLS_CONCAT.out.vcf
        )

        BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT.out.tbi)
            .map { meta, vcf , tbi ->
                def new_meta = meta + [caller:'clairs']
                return [new_meta, vcf, tbi]
            }
            .set{clairs_ch}
        // clairs_ch: [meta(+caller:'clairs'), vcf, tbi]  -- merged and sorted ClairS somatic VCF
    }

    // DEEPSOMATIC: somatic variant calling using deep learning T/N model
    if(params.somatic_var_keep.contains('deepsomatic')) {

        // DeepSomatic expects [normal, tumor] order (opposite of input tuple)
        tumor_normal_bams
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    return [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
            }
            .set{ deepsomatic_input }
        // deepsomatic_input: [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]

        //
        // SUBWORKFLOW: DEEPSOMATIC (local)
        // Input:  [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        //         [[:],[]]  -- empty intervals
        //         fasta / fai / [[:],[]]  -- empty GZI
        // Output: .vcf       -- [meta, vcf]
        //         .vcf_index -- [meta, tbi]
        //
        DEEPSOMATIC (
            deepsomatic_input,
            [[:],[]],  // intervals (empty = genome-wide)
            fasta,
            fai,
            [[:],[]]   // GZI (empty if FASTA is uncompressed)
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
    // If both callers requested: run consensus subworkflow; otherwise pass through single-caller output
    if (params.somatic_var_keep.size() > 1) {
        clairs_ch
            .mix(deepsomatic_ch)
            .set{combine_somatic_ch}
        // combine_somatic_ch: [meta(+caller), vcf, tbi]  -- one item per caller per sample

        // SUBWORKFLOW: SOMATIC_CONSENSUS (SMALL_VARIANT_CONSENSUS alias)
        SOMATIC_CONSENSUS(
            combine_somatic_ch,
            fasta,
            fai,
            params.prioritize_caller_somatic,
            params.somatic_var_combine
        )

        SOMATIC_CONSENSUS.out.vcf
            .join(SOMATIC_CONSENSUS.out.tbi)
            .set{ somatic_vcf }
        // somatic_vcf: [meta(+caller from consensus), vcf, tbi]
    }
    else if (params.somatic_var_keep == ['clair']) {
        clairs_ch
            .set{somatic_vcf}
    }
    else if (params.somatic_var_keep == ['deepsomatic']) {
        deepsomatic_ch
            .set{somatic_vcf}
    }

    // Strip 'caller' from meta before emitting
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

    emit:
    somatic_vcf  // [meta, vcf, tbi]  -- final somatic VCF (ClairS, DeepSomatic, or consensus)
}
