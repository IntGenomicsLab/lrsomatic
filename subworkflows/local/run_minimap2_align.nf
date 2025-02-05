//
// Run minimap2 on tumour (and normal) ubams
//

include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'

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
        .view()
        .set { ch_restructured_ubams }
    // Run minimap2 align on each bam file
    MINIMAP2_ALIGN ( 
        ch_restructured_ubams,
        ch_minimap_index,
        true,
        'bai',
        "", 
        "" 
    )
    
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_aligned = MINIMAP2_ALIGN.out.bam
    
    // Restructure back to original format

        
    emit:
        aligned = ch_aligned
        versions = ch_versions 
}
