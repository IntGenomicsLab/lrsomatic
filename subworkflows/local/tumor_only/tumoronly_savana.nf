// IMPORT MODULES
include { SAVANA_TO         } from '../../../modules/nf-core/savana/to/main'
include { SAVANA_PON_FILTER } from '../../../modules/local/savana_pon_filter/main'

workflow TUMORONLY_SAVANA {

    take:
    tumor_input // [meta, tumor_bam, tumor_bai, phased_germline_vcf, phased_germline_tbi]  -- tumor-only, no matched normal
    fasta       // [[:], fasta]
    fai         // [[:], fai]
    contigs     // [[:], contigs]  -- one contig per line, restricts SAVANA to canonical chromosomes
    g1000_vcf   // [[:], '1000g_hg38' | '1000g_t2t']  -- bundled 1000 Genomes population SNP set (genome-specific)
    sv_pon      // [[:], pon.tsv.gz] or [[:], []]  -- Severus germline SV panel of normals; PacBio calls are filtered against it

    main:
    // .first() re-establishes a value channel: fasta/fai are singletons reused for every sample,
    // but .map() on fai strips that property, so plain .combine() would only fire once for the run.
    ch_fasta_fai = fasta.combine(fai.map { _meta, index -> index }).first()
    // ch_fasta_fai: [[:], fasta, fai]

    //
    // MODULE: SAVANA_TO (label: process_high)
    // Chains breakpoint-calling + classification + CNA internally. Tumor-only has no matched
    // germline control, so allele counting uses the bundled 1000g population SNP set
    // (--g1000_vcf) rather than the tumour's own phased VCF, per SAVANA's README.
    // Input:  [meta, tumor_bam, tumor_bai, [], [], []]  -- snp_vcf/allele_counts_het_snps/breakpoints unused
    //         [[:], fasta, fai] / contigs / [[:],[]] blacklist / g1000_vcf
    // Output: .somatic_vcf -- [meta, vcf]  -- classified somatic SV VCF (savana classify, run internally)
    //         .cna         -- [meta, tsv]  -- segmented absolute copy number (savana cna, run internally)
    //
    tumor_input
        .map { meta, tumor_bam, tumor_bai, _phased_vcf, _phased_tbi ->
            def snp_vcf = []
            def allele_counts_het_snps = []
            def breakpoints = []
            return [meta, tumor_bam, tumor_bai, snp_vcf, allele_counts_het_snps, breakpoints]
        }
        .set { savana_to_input }
    // savana_to_input: [meta, tumor_bam, tumor_bai, [], [], []]

    SAVANA_TO (
        savana_to_input,
        ch_fasta_fai,
        contigs,
        [[:], []],  // blacklist (unused)
        g1000_vcf
    )

    //
    // MODULE: SAVANA_PON_FILTER (label: process_single)
    // SAVANA's PacBio classifier only rejects a call on normal-read evidence, so without a normal every germline
    // SV with enough tumour reads comes out somatic (~28k on a PacBio genome, against ~2 with its matched normal).
    // ONT uses SAVANA's tumour-only model and is left alone. Without a PoN for the genome the VCF passes unfiltered.
    // Input:  [meta, classified.somatic.vcf] / [[:], pon.tsv.gz]
    // Output: .vcf -- [meta, vcf]  -- the somatic VCF without events in the panel of normals
    //
    SAVANA_TO.out.somatic_vcf
        .combine(sv_pon)
        .branch { meta, vcf, _pon_meta, pon ->
            filter: meta.platform == 'pb' && pon
                return [meta, vcf]
            keep: true
                return [meta, vcf]
        }
        .set { savana_somatic_vcf }
    // savana_somatic_vcf.filter / .keep: [meta, vcf]

    SAVANA_PON_FILTER (
        savana_somatic_vcf.filter,
        sv_pon
    )

    emit:
    somatic_vcf          = savana_somatic_vcf.keep.mix(SAVANA_PON_FILTER.out.vcf) // [meta, vcf] -- classified somatic SV VCF (PoN-filtered for PacBio)
    cn_calls             = SAVANA_TO.out.cna                   // [meta, tsv]    -- segmented absolute copy number
    fitted_purity_ploidy = SAVANA_TO.out.fitted_purity_ploidy  // [meta, tsv]    -- selected purity/ploidy fit (absent when no fit)
    allele_counts        = SAVANA_TO.out.allele_counts         // [meta, bed]    -- het-SNP allele counts (ReConPlot BAF track)
}
