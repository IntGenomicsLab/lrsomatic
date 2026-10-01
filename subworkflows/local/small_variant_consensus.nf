include { BCFTOOLS_NORM                                      } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_REJOIN              } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC                                      } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_QUERY                                     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE                                  } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_SORT as SORT_POST_NORM                    } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CONSENSUS           } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_CALLER_UNION                              } from '../../modules/local/bcftools/callerunion/main'



workflow SMALL_VARIANT_CONSENSUS {
    take:
    mixed_vcfs       // [meta(+caller field), vcf, tbi]  -- one item per caller per sample
    //                    meta.caller is one of: 'clair3', 'clairs-to', 'clairs', 'deepvariant', 'deepsomatic'
    fasta            // [[:], fasta]
    _fai             // [[:], fai]
    prioritize_caller // str: which caller's calls take priority ('deepvariant'/'deepsomatic' or 'clair')
    combine_method   // str: 'consensus' (shared calls only) or 'all' (union of both callers' calls)

    main:

    if (!(combine_method in ['consensus', 'all'])) {
        error("combine_method must be 'consensus' or 'all', got '${combine_method}'")
    }
    if (!(prioritize_caller in ['deepvariant', 'deepsomatic', 'clair'])) {
        error("prioritize_caller must be one of [deepvariant, deepsomatic, clair], got '${prioritize_caller}'")
    }

    //
    // MODULE: BCFTOOLS_NORM (label: process_medium) -- left-align; in 'consensus' mode also split multi-allelics for isec
    // Input:  [meta(+split), vcf, tbi]  -- per-caller VCF
    // Output: .vcf -- [meta, vcf]  -- left-aligned, normalised VCF (unsorted)
    //
    BCFTOOLS_NORM(
        mixed_vcfs.map { meta, vcf, tbi -> [meta + [split: combine_method == 'consensus'], vcf, tbi] },
        fasta
    )

    //
    // MODULE: SORT_POST_NORM (BCFTOOLS_SORT alias, label: process_medium) -- re-sort and index after normalisation
    // Input:  [meta, vcf]
    // Output: .vcf -- [meta, vcf.gz]
    //         .tbi -- [meta, tbi]
    //
    SORT_POST_NORM(BCFTOOLS_NORM.out.vcf)

    SORT_POST_NORM.out.vcf
        .join(SORT_POST_NORM.out.tbi)
        .set { normalized_vcfs }
    // normalized_vcfs: [meta(+caller), vcf.gz, tbi]  -- normalised, sorted per-caller VCF

    //
    // ALLELE FREQUENCY KEY -- BCFTOOLS_ANNOTATE below renames the AF FORMAT field to the priority caller's:
    //   FORMAT/AF  -> FORMAT/VAF  when prioritize_caller is 'deepvariant'/'deepsomatic'
    //   FORMAT/VAF -> FORMAT/AF   when prioritize_caller is 'clair'
    // Only 'all' mode renames: it merges both callers, so the merged VCF needs one AF key for WAKHAN.

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
        .join(BCFTOOLS_QUERY.out.output, failOnMismatch: true, failOnDuplicate: true)
        .join(BCFTOOLS_QUERY.out.index, failOnMismatch: true, failOnDuplicate: true)
        .map{ meta, vcf, tbi, annotations, annotations_index ->
                    def columns = []       // no extra column specs
                    def header_lines = []  // no extra header lines
                    def rename_chrs = []   // no chromosome renaming
                    // 'all' mode merges both callers, so unify the AF key; 'consensus' needs no rename.
                    def new_meta = combine_method == 'all'
                        ? meta + [rename_to: (prioritize_caller in ['deepvariant', 'deepsomatic'] ? 'VAF' : 'AF')]
                        : meta
                return [ new_meta, vcf, tbi, annotations, annotations_index, columns, header_lines, rename_chrs ]
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
        .join(BCFTOOLS_ANNOTATE.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi ->
            def clean_meta = meta.findAll { k, _v -> !(k in ['rename_to', 'split']) }
            return [clean_meta, vcf, tbi]
        }
        .set{annotated_vcfs}
    // annotated_vcfs: [meta(+caller), vcf, tbi]  -- VCF with CALLER INFO tag

    // Branch annotated VCFs by caller family for the intersection step
    // `other` errors on an unrecognised meta.caller instead of silently dropping the sample.
    annotated_vcfs
        .branch { meta, _vcfs, _tbi ->
            deepvariant: meta.caller in [ 'deepvariant', 'deepsomatic' ]
            clair: meta.caller in ['clair3','clairs-to','clairs']
            other: true
        }
        .set{annotated_vcfs_branched}

    annotated_vcfs_branched.other
        .map { meta, _vcfs, _tbi ->
            error("SMALL_VARIANT_CONSENSUS: unrecognised meta.caller '${meta.caller}' for sample '${meta.id}'; expected one of [deepvariant, deepsomatic, clair3, clairs-to, clairs]")
        }
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
    // failOnMismatch: a sample missing one caller would otherwise be dropped silently.
    deepvariant_ch
        .join(clair_ch, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
            def vcfs = [deepvar_vcf, clair_vcf]
            def tbis = [deepvar_tbi, clair_tbi]
            return [ meta, vcfs, tbis]
        }
        .set{mixed_vcfs}
    // mixed_vcfs (re-paired): [meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]

    if (combine_method == 'consensus') {
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
        //
        BCFTOOLS_ISEC(isec_input)

        // Take only the intersection: variants called by BOTH callers, from the prioritized caller's record
        def isec_consensus_vcf = prioritize_caller in ['deepvariant', 'deepsomatic']
            ? BCFTOOLS_ISEC.out.deepvar_consensus_vcf
            : BCFTOOLS_ISEC.out.clair_consensus_vcf
        // ISEC always writes 0002.vcf.gz, so germline and somatic would collide by basename in BCFTOOLS_CONCAT;
        // BCFTOOLS_SORT_CONSENSUS renames it per sample (conf/modules.config)
        BCFTOOLS_SORT_CONSENSUS(isec_consensus_vcf)

        //
        // MODULE: BCFTOOLS_NORM_REJOIN (BCFTOOLS_NORM alias) -- rejoin split sites (-m +any) so LongPhase and Wakhan see one record per position
        // Input:  [meta, vcf, tbi]  -- sorted consensus VCF (one caller's records only)
        // Output: .vcf -- [meta, vcf.gz]
        //         .tbi -- [meta, tbi]
        //
        BCFTOOLS_NORM_REJOIN(
            BCFTOOLS_SORT_CONSENSUS.out.vcf.join(BCFTOOLS_SORT_CONSENSUS.out.tbi, failOnMismatch: true, failOnDuplicate: true),
            fasta
        )
        BCFTOOLS_NORM_REJOIN.out.vcf.set{ vcf }
        BCFTOOLS_NORM_REJOIN.out.tbi.set{ tbi }
        // vcf/tbi: [meta, vcf/tbi]  -- consensus calls from the priority caller, multi-allelics rejoined
    }
    else {
        // Union by locus, one caller's record per position: PASS records first, then the priority caller's.
        // Records are unsplit here, so every output record is one caller's call.
        mixed_vcfs
            .map { meta, vcfs, tbis ->
                def prio = prioritize_caller in ['deepvariant', 'deepsomatic'] ? 0 : 1
                return [meta, vcfs[prio], tbis[prio], vcfs[1 - prio], tbis[1 - prio]]
            }
            .set{ union_input }
        // union_input: [meta, prio_vcf, prio_tbi, other_vcf, other_tbi]

        //
        // MODULE: BCFTOOLS_CALLER_UNION (label: process_single)
        // Input:  [meta, prio_vcf, prio_tbi, other_vcf, other_tbi]
        // Output: .vcf / .tbi -- [meta, vcf.gz] / [meta, tbi]  -- sorted union VCF, one record per position
        //
        BCFTOOLS_CALLER_UNION(union_input)
        BCFTOOLS_CALLER_UNION.out.vcf.set{ vcf }
        BCFTOOLS_CALLER_UNION.out.tbi.set{ tbi }
    }

    emit:
    vcf  // [meta, vcf]  -- final consensus/combined VCF, one record per position
    tbi  // [meta, tbi]

}
