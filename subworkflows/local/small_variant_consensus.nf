
include { BCFTOOLS_NORM        } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC        } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_MERGE        } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY        } from '../../modules/nf-core/bcftools/merge/main'


workflow SMALL_VARIANT_CONSENSUS.nf{
    take:
    mixed_vcfs // [meta: w caller_info,mixed_vcfs, mixed_indicies]
    fasta
    fai
    var_keep_method

    main:
    //normalize VCFs
    BCFTOOLS_NORM(mixed_vcfs, fasta)

    BCFTOOLS_NORM.out.vcf   
             .join(BCFTOOLS_NORM.tbi)
             .set {normalized_vcfs}

    // create annotation file with caller name
    BCFTOOLS_QUERY(normalized_vcfs, [], [], [])

    normalized_vcfs
        .join(BCFTOOLS_QUERY.out.output)
        .map{ meta, vcf, tbi, annotations ->
                    def annotations_index = []
                    def columns = []
                    def header_lines = []
                    def rename_chrs = []
                return [ meta, vcf, tbi,file,annotations,annotations_index, columns, header_lines, rename_chrs ]
             }
             .set{annotate_input}

    // Annotate vcfs with caller id
    BCFTOOLS_ANNOTATE(annotate_input)

    BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.tbi)
        .set{annotated_vcfs}

    annotated_vcfs
        .branch { meta, _vcfs, _tbi ->
            deepvariant: meta.caller in [ 'deepvariant', 'deepsomatic' ]
            clair: meta.caller in ['clair3','clairs-to','clairs']
        }
        .set{annotated_vcfs_branched}

    clair_ch = annotated_vcfs_branched.clair
    deepvariant = annotated_vcfs_branched.deepvariant
    
    clair_ch.
        map {meta, vcfs, tbi ->
            meta
        }
    if (var_keep_method == 'consensus') {

        BCFTOOLS_NORM.out.vcf   
             .join(BCFTOOLS_NORM.tbi)
             .map{ meta, vcf, tbi ->
                    def file = []
                    def target = []
                    def regions = []
                return [meta, vcf, tbi, file, target, regions]
             }
             .set{isec_input}
        BCFTOOLS_ISEC(isec_input)
    }
    else if (var_keep_method == 'all'){
        
        BCFTOOLS_ANNOTATE.out.vcf
            .join(BCFTOOLS_ANNOTATE.out.tbi)
            .map{ meta, vcf, tbi ->
                def bed = []
                return [ bed ]
            }
            .set{ merge input}
        fasta
            .join(fai)
            .set{ fasta_fai }
        BCFTOOLS_MERGE(merge_input, fasta_fai)
    }

}