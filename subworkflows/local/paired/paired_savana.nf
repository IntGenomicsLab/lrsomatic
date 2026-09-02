// IMPORT MODULES
include { SAVANA_RUN      } from '../../../modules/nf-core/savana/run/main'
include { SAVANA_CLASSIFY } from '../../../modules/nf-core/savana/classify/main'
include { SAVANA_CNA      } from '../../../modules/nf-core/savana/cna/main'

workflow PAIRED_SAVANA {

    take:
    tumor_normal_input // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, phased_germline_tbi]
    fasta              // [[:], fasta]
    fai                // [[:], fai]

    main:
    ch_fasta_fai = fasta.combine(fai.map { it[1] })
    // ch_fasta_fai: [[:], fasta, fai]

    //
    // MODULE: SAVANA_RUN (label: process_high)
    // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    //         [[:], fasta, fai]
    // Output: .sv_breakpoints_vcf -- [meta, vcf]  -- raw (unclassified) SV breakpoints
    //
    tumor_normal_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _phased_vcf, _phased_tbi ->
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
        .set { savana_run_input }
    // savana_run_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    SAVANA_RUN (
        savana_run_input,
        ch_fasta_fai
    )

    //
    // MODULE: SAVANA_CLASSIFY (label: process_medium)
    // Input:  [meta, sv_breakpoints_vcf]  -- SAVANA_RUN's raw breakpoints
    // Output: .somatic_vcf -- [meta, vcf]  -- classified somatic SV VCF
    //
    SAVANA_CLASSIFY (
        SAVANA_RUN.out.sv_breakpoints_vcf
    )

    //
    // MODULE: SAVANA_CNA (label: process_high)
    // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, snp_vcf, [], breakpoints]
    //           snp_vcf     -- reuses the phased germline VCF already computed for SEVERUS/phasing
    //           breakpoints -- SAVANA_RUN's raw breakpoints (informs CNA binning around SV junctions)
    //         [[:], fasta, fai]
    //         [] / [] / []  -- contigs / blacklist / g1000_vcf (all unused; snp_vcf takes priority)
    // Output: .cna -- [meta, tsv]  -- segmented absolute copy number. Same column shape as the
    //                                 SAVANA-native TSV this pipeline's downstream consumer already
    //                                 parses (chromosome,start,end,segment_id,bin_count,
    //                                 sum_of_bin_lengths,weight,copyNumber[,...] -- columns 1,2,3,8).
    //
    tumor_normal_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, _phased_tbi ->
            def allele_counts_het_snps = []
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, allele_counts_het_snps]
        }
        .join(SAVANA_RUN.out.sv_breakpoints_vcf)
        .set { savana_cna_input }
    // savana_cna_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, snp_vcf, [], breakpoints]

    SAVANA_CNA (
        savana_cna_input,
        ch_fasta_fai,
        [],  // contigs (unused)
        [],  // blacklist (unused)
        []   // g1000_vcf (unused -- snp_vcf takes priority)
    )

    emit:
    somatic_vcf = SAVANA_CLASSIFY.out.somatic_vcf  // [meta, vcf]  -- classified somatic SV VCF
    cn_calls    = SAVANA_CNA.out.cna               // [meta, tsv]  -- segmented absolute copy number
}
