
include { BCFTOOLS_NORM                                      } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC                                      } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_QUERY                                     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE                                  } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as STANDARDIZE_AF                } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_CONCAT                                    } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                      } from '../../modules/nf-core/bcftools/sort/main'



workflow SMALL_VARIANT_CONSENSUS {
    take:
    mixed_vcfs       // [meta(+caller field), vcf, tbi]  -- one item per caller per sample
    //                    meta.caller is one of: 'clair3', 'clairs-to', 'clairs', 'deepvariant', 'deepsomatic'
    fasta            // [[:], fasta]
    _fai             // [[:], fai]
    prioritize_caller // str: which caller's calls take priority ('deepvariant'/'deepsomatic' or 'clair')
    combine_method   // str: 'consensus' (intersection only) or 'all' (intersection + private calls from priority caller)

    main:

    //
    // MODULE: BCFTOOLS_NORM (label: process_medium)
    // Input:  [meta, vcf, tbi]  -- per-caller VCF
    // Output: .vcf -- [meta, vcf]  -- left-aligned, normalised VCF
    //         .tbi -- [meta, tbi]
    //
    BCFTOOLS_NORM(mixed_vcfs, fasta)

    BCFTOOLS_NORM.out.vcf
             .join(BCFTOOLS_NORM.out.tbi)
             .set {normalized_vcfs}
    // normalized_vcfs: [meta(+caller), vcf, tbi]  -- normalised per-caller VCF

    //
    // MODULE: STANDARDIZE_AF (BCFTOOLS_ANNOTATE alias, label: process_low)
    // Renames the allele frequency FORMAT field to match the priority caller's convention:
    //   FORMAT/AF  -> FORMAT/VAF  when prioritize_caller is 'deepvariant'/'deepsomatic'
    //   FORMAT/VAF -> FORMAT/AF   when prioritize_caller is 'clair'
    // This is a no-op for VCFs that already use the target field name.
    //
    normalized_vcfs
        .map { meta, vcf, tbi ->
            def rename_to = prioritize_caller in ['deepvariant', 'deepsomatic'] ? 'VAF' : 'AF'
            def new_meta = meta + [rename_to: rename_to]
            return [new_meta, vcf, tbi, [], [], [], [], []]
        }
        .set { standardize_input }

    STANDARDIZE_AF(standardize_input)

    STANDARDIZE_AF.out.vcf
        .join(STANDARDIZE_AF.out.tbi)
        .map { meta, vcf, tbi ->
            def clean_meta = meta.findAll { k, v -> k != 'rename_to' }
            return [clean_meta, vcf, tbi]
        }
        .set { normalized_vcfs }
    // normalized_vcfs: [meta(+caller), vcf, tbi]  -- normalised, AF-standardized per-caller VCF

    //
    // MODULE: BCFTOOLS_QUERY (label: process_single)
    // Extract variant positions to build a caller-annotation file used by BCFTOOLS_ANNOTATE
    // Input:  [meta, vcf, tbi]  -- normalised VCF
    // Output: .output -- [meta, tsv]  -- tab-separated annotation file (CHROM POS CALLER)
    //         .index  -- [meta, tbi]
    //
    BCFTOOLS_QUERY(normalized_vcfs, [], [], [])

    // Prepare BCFTOOLS_ANNOTATE input: VCF + caller-name annotation file
    normalized_vcfs
        .join(BCFTOOLS_QUERY.out.output)
        .join(BCFTOOLS_QUERY.out.index)
        .map{ meta, vcf, tbi, annotations, annotations_index ->
                    def columns = []       // no extra column specs
                    def header_lines = []  // no extra header lines
                    def rename_chrs = []   // no chromosome renaming
                return [ meta, vcf, tbi, annotations, annotations_index, columns, header_lines, rename_chrs ]
             }
             .set{annotate_input}
    // annotate_input: [meta, vcf, tbi, annotations_tsv, annotations_tbi, [], [], []]

    //
    // MODULE: BCFTOOLS_ANNOTATE (label: process_medium)
    // Adds CALLER INFO field to each VCF record using the query-generated annotation file
    // Input:  [meta, vcf, tbi, annotations_tsv, annotations_tbi, [], [], []]
    // Output: .vcf -- [meta, vcf]  -- VCF with CALLER annotation added
    //         .tbi -- [meta, tbi]
    //
    BCFTOOLS_ANNOTATE(annotate_input)

    BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)
        .set{annotated_vcfs}
    // annotated_vcfs: [meta(+caller), vcf, tbi]  -- VCF with CALLER INFO tag

    // Branch annotated VCFs by caller family for the intersection step
    annotated_vcfs
        .branch { meta, _vcfs, _tbi ->
            deepvariant: meta.caller in [ 'deepvariant', 'deepsomatic' ]
            clair: meta.caller in ['clair3','clairs-to','clairs']
        }
        .set{annotated_vcfs_branched}
    // annotated_vcfs_branched.deepvariant: [meta(caller=deepvariant/deepsomatic), vcf, tbi]
    // annotated_vcfs_branched.clair:       [meta(caller=clair3/clairs-to/clairs), vcf, tbi]

    clair_ch = annotated_vcfs_branched.clair
    deepvariant_ch = annotated_vcfs_branched.deepvariant

    // Strip 'caller' field from meta before joining so both channels share the same key
    clair_ch
        .map {meta, vcfs, tbi ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'type',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return [ new_meta, vcfs, tbi]
        }
        .set{clair_ch}
    // clair_ch: [meta (no caller), vcf, tbi]

    deepvariant_ch
        .map {meta, vcfs, tbi ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'type',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return [ new_meta, vcfs, tbi]
        }
        .set{deepvariant_ch}
    // deepvariant_ch: [meta (no caller), vcf, tbi]

    // Join DeepVariant and Clair VCFs per sample into a single tuple for BCFTOOLS_ISEC
    deepvariant_ch
        .join(clair_ch)
        .map { meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
            def vcfs = [deepvar_vcf, clair_vcf]
            def tbis = [deepvar_tbi, clair_tbi]
            return [ meta, vcfs, tbis]
        }
        .set{mixed_vcfs}
    // mixed_vcfs (re-paired): [meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]

    // Add empty optional fields required by BCFTOOLS_ISEC
    mixed_vcfs
         .map{ meta, vcfs, tbis ->
                def file = []    // no regions file
                def target = []  // no target sites
                def regions = [] // no region string
            return [meta, vcfs, tbis, file, target, regions]
         }
         .set{isec_input}
    // isec_input: [meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi], [], [], []]

    //
    // MODULE: BCFTOOLS_ISEC (label: process_medium)
    // Computes the intersection and private sets for the two callers
    // Input:  [meta, [vcf1, vcf2], [tbi1, tbi2], [], [], []]
    // Output (custom nf-core module outputs):
    //   .deepvar_consensus_vcf  -- [meta, vcf]  -- variants called by both callers (DeepVariant record)
    //   .clair_consensus_vcf    -- [meta, vcf]  -- variants called by both callers (Clair record)
    //   .deepvar_private_vcf    -- [meta, vcf]  -- variants unique to DeepVariant
    //   .clair_private_vcf      -- [meta, vcf]  -- variants unique to Clair
    //   (+ corresponding .tbi outputs for each)
    //
    BCFTOOLS_ISEC(isec_input)

    if (combine_method == 'consensus') {
        // Take only the intersection: variants called by BOTH callers
        // Use the record from the prioritized caller
        if (prioritize_caller in ['deepvariant', 'deepsomatic']) {
            BCFTOOLS_ISEC.out.deepvar_consensus_vcf
                .set{vcf}
            BCFTOOLS_ISEC.out.deepvar_consensus_tbi
                .set{tbi}
        }
        else if (prioritize_caller == 'clair') {
            BCFTOOLS_ISEC.out.clair_consensus_vcf
                .set{vcf}
            BCFTOOLS_ISEC.out.clair_consensus_tbi
                .set{tbi}
        }
        // vcf/tbi: [meta, vcf/tbi]  -- consensus-only calls from the priority caller
    }

    else if (combine_method == 'all') {
        // Take the intersection PLUS the private calls from the prioritized caller
        // (private calls from the non-priority caller are discarded)
        if (prioritize_caller in ['deepvariant', 'deepsomatic']) {
            // consensus (DeepVariant record) + DeepVariant-private variants
            BCFTOOLS_ISEC.out.deepvar_consensus_vcf
                .join(BCFTOOLS_ISEC.out.deepvar_consensus_tbi)
                .join(BCFTOOLS_ISEC.out.clair_private_vcf)
                .join(BCFTOOLS_ISEC.out.clair_private_tbi)
                .map{ meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
            // concat_input: [meta, [consensus_vcf, private_vcf], [consensus_tbi, private_tbi]]
            BCFTOOLS_CONCAT(concat_input)
            BCFTOOLS_CONCAT.out.vcf
                .set{concat_out}
        }
        else if (prioritize_caller == 'clair') {
            // consensus (Clair record) + Clair-private variants
            BCFTOOLS_ISEC.out.deepvar_private_vcf
                .join(BCFTOOLS_ISEC.out.deepvar_private_tbi)
                .join(BCFTOOLS_ISEC.out.clair_consensus_vcf)
                .join(BCFTOOLS_ISEC.out.clair_consensus_tbi)
                .map{ meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
            // concat_input: [meta, [private_vcf, consensus_vcf], [private_tbi, consensus_tbi]]
            BCFTOOLS_CONCAT(concat_input)
            BCFTOOLS_CONCAT.out.vcf
                .set{concat_out}
        }
        // concat_out: [meta, vcf]  -- unsorted concatenated VCF (consensus + priority-caller-private)
        BCFTOOLS_SORT(concat_out)
        BCFTOOLS_SORT.out.vcf
            .set{vcf}
        BCFTOOLS_SORT.out.tbi
            .set{tbi}
        // vcf/tbi: [meta, vcf/tbi]  -- sorted combined VCF
    }

    emit:
    vcf  // [meta, vcf]  -- final consensus/combined VCF
    tbi  // [meta, tbi]

}
