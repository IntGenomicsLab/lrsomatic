//
// Run minimap2 on tumour (and normal) ubams
//

include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TUMOUR } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_NORMAL } from '../modules/nf-core/minimap2/align/main'

workflow RUN_MINIMAP2 {
    take:
        ch_cat_ubams

    main:
        ch_versions = Channel.empty()
				
		// Check if sample is paired or not
		if ch_cat_ubams.meta.paired {
			ubam_normal = ch_cat_ubams.map{ meta, tumor_bam, normal_bam -> [ meta, normal_bam ] }
			ubam_tumor  = ch_cat_ubams.map{ meta, tumor_bam, normal_bam -> [ meta, tumor_bam ] }

			
			MINIMAP2_ALIGN_TUMOUR( ubam_normal )
				.reads
				.set { ch_align_tumour }
				
			MINIMAP2_ALIGN_NORMAL( ubam_tumor )
				.reads
				.set { ch_align_normal }
				
			// Remake channel
			ch_aligned = ch_align_tumour.merge( ch_align_normal )
			//TODO: check if this works
		} else {
			MINIMAP2_ALIGN_TUMOUR( ch_cat_ubams )
				.reads
				.set { ch_aligned }
		}
}