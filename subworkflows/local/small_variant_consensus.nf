include { BCFTOOLS_NORM                                      } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC                                      } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_QUERY                                     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE                                  } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as STANDARDIZE_AF                } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_CONCAT                                    } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                      } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as SORT_POST_NORM                    } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CONSENSUS           } from '../../modules/nf-core/bcftools/sort/main'



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
    // MODULE: BCFTOOLS_NORM (label: process_medium) -- left-align and normalise; sorted after, since left-alignment can reorder records
    // Input:  [meta, vcf, tbi]  -- per-caller VCF
    // Output: .vcf -- [meta, vcf]  -- left-aligned, normalised VCF (unsorted)
    //
    BCFTOOLS_NORM(mixed_vcfs, fasta)

    //
    // MODULE: SORT_POST_NORM (BCFTOOLS_SORT alias, label: process_medium) -- re-sort and index after normalisation
    // Input:  [meta, vcf]
    // Output: .vcf -- [meta, vcf.gz]
    //         .tbi -- [meta, tbi]
    //
    SORT_POST_NORM(BCFTOOLS_NORM.out.vcf)

    SORT_POST_NORM.out.vcf
        .join(SORT_POST_NORM.out.index)
        .set { normalized_vcfs }
    // normalized_vcfs: [meta(+caller), vcf.gz, tbi]  -- normalised, sorted per-caller VCF

    //
    // MODULE: STANDARDIZE_AF (BCFTOOLS_ANNOTATE alias, label: process_low) -- rename the AF FORMAT field to the priority caller's:
    //   FORMAT/AF  -> FORMAT/VAF  when prioritize_caller is 'deepvariant'/'deepsomatic'
    //   FORMAT/VAF -> FORMAT/AF   when prioritize_caller is 'clair'
    //   Only the other family's VCFs are renamed: Clair callers write AF, DeepVariant/DeepSomatic write VAF,
    //   and bcftools >= 1.24 fails on a rename whose source tag is absent instead of skipping it
    //
    if (combine_method == 'all') {
        def rename_to = prioritize_caller in ['deepvariant', 'deepsomatic'] ? 'VAF' : 'AF'
        normalized_vcfs
            .branch { meta, _vcf, _tbi ->
                def caller_tag = meta.caller in ['deepvariant', 'deepsomatic'] ? 'VAF' : 'AF'
                rename: caller_tag != rename_to
                keep: true
            }
            .set { af_tag }

        af_tag.rename
            .map { meta, vcf, tbi ->
                def new_meta = meta + [rename_to: rename_to]
                return [new_meta, vcf, tbi, [], [], [], [], []]
            }
            .set { standardize_input }

        STANDARDIZE_AF(standardize_input)

        STANDARDIZE_AF.out.vcf
            .join(STANDARDIZE_AF.out.index)
            .map { meta, vcf, tbi ->
                def clean_meta = meta.findAll { k, _v -> k != 'rename_to' }
                return [clean_meta, vcf, tbi]
            }
            .mix(af_tag.keep)
            .set { normalized_vcfs }
        // normalized_vcfs: [meta(+caller), vcf, tbi]  -- normalised, AF-standardized per-caller VCF
    }
    // In 'consensus' mode, normalized_vcfs comes from SORT_POST_NORM (post-BCFTOOLS_NORM re-sorting)

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
        .join(BCFTOOLS_ANNOTATE.out.index)
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
    // MODULE: BCFTOOLS_ISEC (label: process_medium) -- shared and private sets of the two callers
    // Input:  [meta, [vcf1, vcf2], [tbi1, tbi2], [], [], []]
    // Output: .deepvar_consensus_vcf / .clair_consensus_vcf -- [meta, vcf]  -- shared calls, DeepVariant or Clair record
    //         .deepvar_private_vcf / .clair_private_vcf     -- [meta, vcf]  -- caller-private calls (+ .tbi for each)
    //
    BCFTOOLS_ISEC(isec_input)

    if (combine_method == 'consensus') {
        // Take only the intersection: variants called by BOTH callers
        // Use the record from the prioritized caller
        if (prioritize_caller in ['deepvariant', 'deepsomatic']) {
            BCFTOOLS_ISEC.out.deepvar_consensus_vcf
                .set{isec_consensus_vcf}
        }
        else if (prioritize_caller == 'clair') {
            BCFTOOLS_ISEC.out.clair_consensus_vcf
                .set{isec_consensus_vcf}
        }
        // ISEC always writes 0002.vcf.gz, so germline and somatic would collide by basename in BCFTOOLS_CONCAT;
        // BCFTOOLS_SORT_CONSENSUS renames it per sample (conf/modules.config)
        BCFTOOLS_SORT_CONSENSUS(isec_consensus_vcf)
        BCFTOOLS_SORT_CONSENSUS.out.vcf.set{vcf}
        BCFTOOLS_SORT_CONSENSUS.out.index.set{tbi}
        // vcf/tbi: [meta, vcf/tbi]  -- consensus-only calls from the priority caller, renamed
    }

    else if (combine_method == 'all') {
        // Intersection plus the priority caller's private calls
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
        BCFTOOLS_SORT.out.index
            .set{tbi}
        // vcf/tbi: [meta, vcf/tbi]  -- sorted combined VCF
    }

    emit:
    vcf  // [meta, vcf]  -- final consensus/combined VCF
    tbi  // [meta, tbi]

}
