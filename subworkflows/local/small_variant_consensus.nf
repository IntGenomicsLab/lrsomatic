include { BCFTOOLS_NORM                                      } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_REJOIN              } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC                                      } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_QUERY                                     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE                                  } from '../../modules/nf-core/bcftools/annotate/main'
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
    combine_method   // str: 'consensus' (shared calls only) or 'all' (union of both callers' calls)

    main:

    //
    // MODULE: BCFTOOLS_NORM (label: process_medium) -- left-align and split multi-allelics for isec; rejoined before phasing
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
        .join(SORT_POST_NORM.out.tbi)
        .set { normalized_vcfs }
    // normalized_vcfs: [meta(+caller), vcf.gz, tbi]  -- normalised, sorted per-caller VCF

    //
    // ALLELE FREQUENCY KEY -- BCFTOOLS_ANNOTATE below renames the AF FORMAT field to the priority caller's:
    //   FORMAT/AF  -> FORMAT/VAF  when prioritize_caller is 'deepvariant'/'deepsomatic'
    //   FORMAT/VAF -> FORMAT/AF   when prioritize_caller is 'clair'
    // This guarantees the merged VCF exposes allele frequency under a single FORMAT key, which is
    // what WAKHAN consumes. Every caller currently emits FORMAT/AF and none emits VAF (verified
    // against Clair3, ClairS-TO, DeepVariant and DeepSomatic output), so under the default
    // prioritize_caller='clair' it is a no-op; it is kept as the guarantee, not the mechanism.
    //
    // The callers do disagree on the AF *declaration*: ClairS-TO says Number=1, DeepVariant and
    // DeepSomatic say Number=A. bcftools concat only warns and keeps the first file's definition.
    // Renaming cannot fix that and `annotate -h` cannot override an existing FORMAT definition,
    // but it is harmless because BCFTOOLS_NORM now splits multi-allelics (-m -any): every record
    // reaching here carries one ALT and one AF value, making the two declarations equivalent.
    // After BCFTOOLS_NORM_REJOIN, consensus files keep the source caller's own header (lossless round
    // trip); 'all' mode takes the Number=A header of isec's 0000.vcf.gz, which concat lists first.
    //
    // The rename is carried out by BCFTOOLS_ANNOTATE below rather than by a second annotate call:
    // meta.rename_to selects the --rename-annots file in conf/modules.config. Only 'all' mode needs
    // it, since in 'consensus' mode every surviving record comes from one caller. BCFTOOLS_QUERY
    // reads only CHROM/POS/REF/ALT, so it does not care whether the rename has happened yet.

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
                    // 'all' mode merges records from both callers into one VCF, so the allele
                    // frequency key is unified here; 'consensus' mode needs no rename.
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
            def clean_meta = meta.findAll { k, _v -> k != 'rename_to' }
            return [clean_meta, vcf, tbi]
        }
        .set{annotated_vcfs}
    // annotated_vcfs: [meta(+caller), vcf, tbi]  -- VCF with CALLER INFO tag

    // Branch annotated VCFs by caller family for the intersection step
    // An unrecognised meta.caller would silently vanish without the `other` arm, taking the whole
    // sample out of the results with a successful exit.
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
    // failOnMismatch: a sample present for one caller but not the other would otherwise be dropped
    // from every downstream result while the run still reported success.
    deepvariant_ch
        .join(clair_ch, failOnMismatch: true, failOnDuplicate: true)
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
        else {
            error("prioritize_caller must be one of [deepvariant, deepsomatic, clair], got '${prioritize_caller}'")
        }
        // ISEC always writes 0002.vcf.gz, so germline and somatic would collide by basename in BCFTOOLS_CONCAT;
        // BCFTOOLS_SORT_CONSENSUS renames it per sample (conf/modules.config)
        BCFTOOLS_SORT_CONSENSUS(isec_consensus_vcf)
        BCFTOOLS_SORT_CONSENSUS.out.vcf.set{vcf}
        BCFTOOLS_SORT_CONSENSUS.out.tbi.set{tbi}
        // vcf/tbi: [meta, vcf/tbi]  -- consensus-only calls from the priority caller, renamed
    }

    else if (combine_method == 'all') {
        // Union: all variants from both callers. Shared variants contribute a single record, taken
        // from the prioritized caller; both callers' private calls are kept. prioritize_caller only
        // selects whose record is used for shared variants, never which calls are kept.
        // The three isec sets are disjoint by construction, so BCFTOOLS_CONCAT needs no -d.
        if (prioritize_caller in ['deepvariant', 'deepsomatic']) {
            // shared (DeepVariant record) + DeepVariant-private + Clair-private
            BCFTOOLS_ISEC.out.deepvar_consensus_vcf
                .join(BCFTOOLS_ISEC.out.deepvar_consensus_tbi)
                .join(BCFTOOLS_ISEC.out.deepvar_private_vcf)
                .join(BCFTOOLS_ISEC.out.deepvar_private_tbi)
                .join(BCFTOOLS_ISEC.out.clair_private_vcf)
                .join(BCFTOOLS_ISEC.out.clair_private_tbi)
                .map{ meta, shared_vcf, shared_tbi, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [shared_vcf, deepvar_vcf, clair_vcf], [shared_tbi, deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
        }
        else if (prioritize_caller == 'clair') {
            // shared (Clair record) + DeepVariant-private + Clair-private
            BCFTOOLS_ISEC.out.clair_consensus_vcf
                .join(BCFTOOLS_ISEC.out.clair_consensus_tbi)
                .join(BCFTOOLS_ISEC.out.deepvar_private_vcf)
                .join(BCFTOOLS_ISEC.out.deepvar_private_tbi)
                .join(BCFTOOLS_ISEC.out.clair_private_vcf)
                .join(BCFTOOLS_ISEC.out.clair_private_tbi)
                .map{ meta, shared_vcf, shared_tbi, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [shared_vcf, deepvar_vcf, clair_vcf], [shared_tbi, deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
        }
        else {
            error("prioritize_caller must be one of [deepvariant, deepsomatic, clair], got '${prioritize_caller}'")
        }
        // concat_input: [meta, [shared_vcf, deepvar_private_vcf, clair_private_vcf], [tbis...]]
        BCFTOOLS_CONCAT(concat_input)
        BCFTOOLS_CONCAT.out.vcf
            .set{concat_out}
        // concat_out: [meta, vcf]  -- unsorted union of both callers' calls
        BCFTOOLS_SORT(concat_out)
        BCFTOOLS_SORT.out.vcf
            .set{vcf}
        BCFTOOLS_SORT.out.tbi
            .set{tbi}
        // vcf/tbi: [meta, vcf/tbi]  -- sorted union VCF
    }

    else {
        error("combine_method must be 'consensus' or 'all', got '${combine_method}'")
    }

    //
    // MODULE: BCFTOOLS_NORM_REJOIN (BCFTOOLS_NORM alias) -- rejoin split sites (-m +any) so LongPhase and Wakhan see one record per position
    // Input:  [meta, vcf, tbi]  -- sorted consensus/union VCF
    // Output: .vcf -- [meta, vcf.gz]
    //         .tbi -- [meta, tbi]
    //
    BCFTOOLS_NORM_REJOIN(
        vcf.join(tbi, failOnMismatch: true, failOnDuplicate: true),
        fasta
    )

    BCFTOOLS_NORM_REJOIN.out.vcf
        .join(BCFTOOLS_NORM_REJOIN.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, rejoined_vcf, rejoined_tbi ->
            vcf: [meta, rejoined_vcf]
            tbi: [meta, rejoined_tbi]
        }
        .set { rejoined }

    emit:
    vcf = rejoined.vcf  // [meta, vcf]  -- final consensus/combined VCF, multi-allelics rejoined
    tbi = rejoined.tbi  // [meta, tbi]

}
