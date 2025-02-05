//
// Run minimap2 on tumour (and normal) ubams
//

include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PB  } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ONT } from '../../modules/nf-core/minimap2/align/main'

workflow RUN_MINIMAP2_ALIGN {
    take:
    ch_ubams // channel: val(meta), [ bam_tumor, bam_normal ] or val(meta), [ bam_tumor ]
    ch_minimap_index

            
    main:
    ch_versions = Channel.empty()
    ch_aligned = Channel.empty()  // Initialize ch_aligned
    
    // Extract all bam files from the channel, keeping it linked to metadata
    ch_ubams
        .flatMap { meta, tumor_bam, normal_bam ->
            def meta_tumor = meta.clone()
            meta_tumor.type = 'tumor'
            def result = [[meta_tumor, tumor_bam]]
            
            if (normal_bam) {
                def meta_normal = meta.clone()
                meta_normal.type = 'normal'
                result << [meta_normal, normal_bam]
            }
            
            return result
        }
        .set { ch_restructured_ubams }
        
    // Split the channel into pacbio and ont
    ch_split = ch_restructured_ubams
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
    
    // Mix back in together
    ch_aligned = ch_aligned_pb.mix(ch_aligned_ont)
    
    // TODO: Restructure back to normal if needed
    // ch_aligned is now [[meta], [bam]] 
    // With meta consisting of [id, paired_data, method, specs, type]
    
    emit:
        aligned = ch_aligned
        versions = ch_versions 
}
