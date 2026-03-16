
include { BCFTOOLS_NORM            } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC            } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_MERGE           } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY           } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_ANNOTATE        } from '../../modules/nf-core/bcftools/annotate/main'



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

    deepvariant_ch
        .join(clair_ch)
        .map { meta, deepvar_vcf, deepvar_tbi, clair_vcf, clair_tbi ->
            def vcfs = [deepvar_vcf, clair_vcf]
            def tbis = [deepvar_tbi, clair_tbi]
            return [ meta, vcfs, tbis]
        }
        .set{mixed_vcfs}
    
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
        BCFTOOLS_ISEC.out.deepvar_style_consensus_vcf
            .set{vcf}
        BCFTOOLS_ISEC.out.deepvar_style_consensus_tbi
            .set{tbi}
    }

    else if (var_keep_method == 'all'){
        mixed_vcfs
            .map{ meta, vcfs, tbis ->
                def bed = []
                return [ meta, vcfs, tbis, bed ]
            }
            .set{ merge_input}
        fasta
            .join(fai)
            .set{ fasta_fai }
        BCFTOOLS_MERGE(merge_input, fasta_fai)
        BCFTOOLS_MERGE.out.vcf
            .set{vcf}
        BCFTOOLS_MERGE.out.index
            .set{tbi}
    }

    emit:
    vcf
    tbi

}