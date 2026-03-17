
include { BCFTOOLS_NORM                                      } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC                                      } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_QUERY                                     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE                                  } from '../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_CONCAT                                    } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                      } from '../../modules/nf-core/bcftools/sort/main'



workflow SMALL_VARIANT_CONSENSUS {
    take:
    mixed_vcfs // [meta: w caller_info,mixed_vcfs, mixed_indicies]
    fasta
    fai
    var_keep_method

    main:
    //normalize VCFs
    BCFTOOLS_NORM(mixed_vcfs, fasta)

    BCFTOOLS_NORM.out.vcf   
             .join(BCFTOOLS_NORM.out.tbi)
             .set {normalized_vcfs}

    // create annotation file with caller name
    BCFTOOLS_QUERY(normalized_vcfs, [], [], [])

    normalized_vcfs
        .join(BCFTOOLS_QUERY.out.output)
        .join(BCFTOOLS_QUERY.out.index)
        .map{ meta, vcf, tbi, annotations, annotations_index ->
                    def columns = []
                    def header_lines = []
                    def rename_chrs = []
                return [ meta, vcf, tbi, annotations, annotations_index, columns, header_lines, rename_chrs ]
             }
             .set{annotate_input}

    // Annotate vcfs with caller id
    BCFTOOLS_ANNOTATE(annotate_input)

    BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)
        .set{annotated_vcfs}
    annotated_vcfs
        .branch { meta, _vcfs, _tbi ->
            deepvariant: meta.caller in [ 'deepvariant', 'deepsomatic' ]
            clair: meta.caller in ['clair3','clairs-to','clairs']
        }
        .set{annotated_vcfs_branched}

    clair_ch = annotated_vcfs_branched.clair
    deepvariant_ch = annotated_vcfs_branched.deepvariant
    
    clair_ch.
        map {meta, vcfs, tbi ->
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
    deepvariant_ch.view()
    clair_ch.view()
    deepvariant_ch
        .join(clair_ch)
        .map { meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
            def vcfs = [deepvar_vcf, clair_vcf]
            def tbis = [deepvar_tbi, clair_tbi]
            return [ meta, vcfs, tbis]
        }
        .set{mixed_vcfs}

    mixed_vcfs.view()
    if (var_keep_method == 'consensus') {
        mixed_vcfs
             .map{ meta, vcfs, tbis ->
                    def file = []
                    def target = []
                    def regions = []
                return [meta, vcfs, tbis, file, target, regions]
             }
             .set{isec_input}
        BCFTOOLS_ISEC(isec_input)

        if (params.trust_caller = 'deepvariant') {
            BCFTOOLS_ISEC.out.clair_consensus_vcf
            .set{vcf}
            BCFTOOLS_ISEC.out.clair_consensus_tbi
            .set{tbi}
        }
        if (params.trust_caller = 'clair') {
            BCFTOOLS_ISEC.out.clair_consensus_vcf
            .set{vcf}
            BCFTOOLS_ISEC.out.clair_consensus_tbi
            .set{tbi}
        }
        
    }

    else if (var_keep_method == 'all'){
        
        mixed_vcfs
             .map{ meta, vcfs, tbis ->
                    def file = []
                    def target = []
                    def regions = []
                return [meta, vcfs, tbis, file, target, regions]
             }
             .set{isec_input}
        
        BCFTOOLS_ISEC(isec_input)

        if (params.trust_caller = 'deepvariant') {
            BCFTOOLS_ISEC.out.deepvar_consensus_vcf
                .join(BCFTOOLS_ISEC.out.deepvar_consensus_tbi)
                .join(BCFTOOLS_ISEC.out.clair_private_vcf)
                .join(BCFTOOLS_ISEC.out.clair_private_tbi)
                .map{ meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
            BCFTOOLS_CONCAT(concat_input)
            BCFTOOLS_CONCAT.out.vcf
                .join(BCFTOOLS_CONCAT.out.tbi)
                .set{concat_out}
        }

        else if (params.trust_caller = 'clair') {
            BCFTOOLS_ISEC.out.deepvar_private_vcf
                .join(BCFTOOLS_ISEC.out.deepvar_private_tbi)
                .join(BCFTOOLS_ISEC.out.clair_consensus_vcf)
                .join(BCFTOOLS_ISEC.out.clair_consensus_tbi)
                .map{ meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
                        return[meta, [deepvar_vcf, clair_vcf], [deepvar_tbi, clair_tbi]]
                }
                .set{concat_input}
            BCFTOOLS_CONCAT(concat_input)
            BCFTOOLS_CONCAT.out.vcf
                .join(BCFTOOLS_CONCAT.tbi)
                .set{concat_out}
        }
        concat_out.view()
        BCFTOOLS_SORT(concat_out)
        BCFTOOLS_SORT.out.vcf
            .set{vcf}
        BCFTOOLS_SORT.out.tbi
            .set{tbi}
    }

    emit:
    vcf
    tbi

}