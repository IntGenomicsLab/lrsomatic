include {ENSEMBLVEP_DOWNLOAD } from '../../modules/nf-core/ensemblvep/download/main.nf'

workflow PREPARE_ANNOTATION {

	take:
		vep_cache
		vep_cache_version
		vep_genome
		vep_args
		vep_species
		download_vep_cache

	main:

		ch_versions = Channel.empty()
		ensemblvep_cache = Channel.empty()
		
		if (download_vep_cache) {
			vep_download_info = Channel.of([[],vep_genome, vep_species, vep_cache_version])
			ENSEMBLVEP_DOWNLOAD(vep_download_info)
			ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
			ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

		}
		else if (!vep_cache) {
			def vep_annotation_cache_key = (vep_cache == "s3://annotation-cache/vep_cache/") ? "${vep_cache_version}_${vep_genome}/" : ""
			def vep_species_suffix = vep_args.contains("--merged") ? '_merged' : (vep_args.contains("--refseq") ? '_refseq' : '')
			def vep_cache_dir = "${vep_annotation_cache_key}${vep_species}${vep_species_suffix}/${vep_cache_version}_${vep_genome}"
			def vep_cache_path_full = file("$vep_cache/$vep_cache_dir", type: 'dir')
			if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
				if (vep_cache == "s3://annotation-cache/vep_cache/") {
					error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
				} else {
					error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}./n")
				}
			}

			ensemblvep_cache = Channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
		}



	//
	// MODULE: ENSEMBLVEP_DOWNLOAD
	//


	emit:
		vep_cache = ensemblvep_cache
		versions = ch_versions 

}