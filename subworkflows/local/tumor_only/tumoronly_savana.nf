// IMPORT MODULES
include { SAVANA_TO } from '../../../modules/nf-core/savana/to/main'

workflow TUMORONLY_SAVANA {

    take:
    tumor_input // [meta, tumor_bam, tumor_bai, phased_germline_vcf, phased_germline_tbi]  -- tumor-only, no matched normal
    fasta       // [[:], fasta]
    fai         // [[:], fai]

    main:
    ch_fasta_fai = fasta.combine(fai.map { it[1] })
    // ch_fasta_fai: [[:], fasta, fai]

    //
    // MODULE: SAVANA_TO (label: process_high)
    // `savana to` chains breakpoint-calling + classification + CNA internally in one call --
    // SAVANA's own docs recommend the matched tumor/normal mode for best performance and flag
    // tumor-only as a fallback; we still wire it in (mirroring how this pipeline already supports
    // tumor-only variants of its other callers) and let users opt in.
    //
    // Input:  [meta, tumor_bam, tumor_bai, snp_vcf, [], []]
    //           snp_vcf -- reuses the phased germline VCF already computed for SEVERUS/phasing; the
    //                      same allele-frequency source PAIRED_SAVANA's SAVANA_CNA call uses
    //         [[:], fasta, fai]
    //         [[:],[]] / [[:],[]] / [[:],[]]  -- contigs / blacklist / g1000_vcf (all unused)
    // Output: .somatic_vcf -- [meta, vcf]  -- classified somatic SV VCF (savana classify, run internally)
    //         .cna         -- [meta, tsv]  -- segmented absolute copy number (savana cna, run internally),
    //                                         optional -- only produced when allele-frequency info
    //                                         (snp_vcf here) is available
    //
    tumor_input
        .map { meta, tumor_bam, tumor_bai, phased_vcf, _phased_tbi ->
            def allele_counts_het_snps = []
            def breakpoints = []
            return [meta, tumor_bam, tumor_bai, phased_vcf, allele_counts_het_snps, breakpoints]
        }
        .set { savana_to_input }
    // savana_to_input: [meta, tumor_bam, tumor_bai, snp_vcf, [], []]

    SAVANA_TO (
        savana_to_input,
        ch_fasta_fai,
        [[:], []],  // contigs (unused)
        [[:], []],  // blacklist (unused)
        [[:], []]   // g1000_vcf (unused -- snp_vcf takes priority)
    )

    emit:
    somatic_vcf = SAVANA_TO.out.somatic_vcf  // [meta, vcf]  -- classified somatic SV VCF
    cn_calls    = SAVANA_TO.out.cna          // [meta, tsv]  -- segmented absolute copy number
}
