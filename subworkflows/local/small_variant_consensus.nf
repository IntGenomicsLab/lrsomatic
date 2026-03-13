
include { BCFTOOLS_NORM        } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ISEC        } from '../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_MERGE        } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY        } from '../../modules/nf-core/bcftools/merge/main'


workflow SMALL_VARIANT_CONSENSUS.nf{
    take:
    ch_vcf // [meta, vcfs, tbis]
    fasta
    fai
    var_keep_method

    main:

    if (var_keep_method == 'consensus') {
        BCFTOOLS_NORM(ch_vcf, fasta)
        BCFTOOLS_NORM.out.vcf   
             .join(BCFTOOLS_NORM.tbi)
             .map{ meta, vcf, tbi ->
                    def file = []
                    def target = []
                    def regions = []
                return [meta, vcf, tbi,file,target,regions]
             }
             .set{isec_input}
        BCFTOOLS_ISEC(isec_input)
    }
    else if (var_keep_method == 'all'){
        BCFTOOLS_NORM(ch_vcf,fasta)
        BCFTOOLS_NORM.out.vcf   
             .join(BCFTOOLS_NORM.tbi)
             .set { vcf_tbi}
        BCFTOOLS_QUERY(vcf_tbi, [], [], [])
        vcf_tbi
        .join(BCFTOOLS_QUERY.out.output)
        .map{ meta, vcf, tbi, annotations ->
                    def annotations_index = []
                    def columns = []
                    def header_lines = []
                    def rename_chrs = []
                return [ meta, vcf, tbi,file,annotations,annotations_index, columns, header_lines, rename_chrs ]
             }
             .set{annotate_input}

        BCFTOOLS_ANNOTATE(ch_vcf)
        
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