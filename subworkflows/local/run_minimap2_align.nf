//
// Run minimap2 on tumour (and normal) ubams
//

include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PB  } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ONT } from '../../modules/nf-core/minimap2/align/main'

workflow RUN_MINIMAP2_ALIGN {
    take:
    ch_ubams // channel: val(meta), [ bam ]
    ch_minimap_index

            
    main:
    ch_versions = Channel.empty()
        
    // Split the channel into pacbio and ont
    ch_split = ch_ubams
        .branch { meta, bam -> 
                ont: meta.method == 'ont'
                pb: meta.method == 'pb'
        }

    // Run minimap2 on PacBio samples
    MINIMAP2_ALIGN_PB ( 
        ch_split.pb,
        ch_minimap_index,
        true,
        'bai',
        "", 
        "" 
    )
    
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_PB.out.versions)
    ch_aligned_pb = MINIMAP2_ALIGN_PB.out.bam
    ch_index_pb = MINIMAP2_ALIGN_PB.out.index


    // Run minimap2 on ONT samples
    MINIMAP2_ALIGN_ONT ( 
        ch_split.ont,
        ch_minimap_index,
        true,
        'bai',
        "", 
        "" 
    )
    
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_ONT.out.versions)
    ch_aligned_ont = MINIMAP2_ALIGN_ONT.out.bam
    ch_index_ont = MINIMAP2_ALIGN_ONT.out.index
    
    // Mix back in together
    ch_aligned = ch_aligned_pb.mix(ch_aligned_ont)
    ch_index = ch_index_pb.mix(ch_index_ont)
    
    emit:
        aligned = ch_aligned
        index = ch_index
        versions = ch_versions 
}
