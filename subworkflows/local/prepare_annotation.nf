include {ENSEMBLVEP_DOWNLOAD } from '../../modules/nf-core/ensemblvep/download/main.nf'

workflow PREPARE_ANNOTATION {

	take:
		vep_cache         // path: local VEP cache directory (or S3 annotation-cache URL)
		vep_cache_version // int:  VEP cache version (e.g. 110)
		vep_genome        // str:  genome assembly string (e.g. "GRCh38")
		vep_args          // str:  extra VEP CLI arguments (parsed to detect --merged / --refseq)
		vep_species       // str:  species name (e.g. "homo_sapiens")
		download_vep_cache // bool: if true, download cache via ENSEMBLVEP_DOWNLOAD instead of using local path

	main:

		ch_versions = channel.empty()
		ensemblvep_cache = channel.empty()

		//
		// MODULE: ENSEMBLVEP_DOWNLOAD (label: process_medium)
		// Only runs when params.download_vep_cache == true
		// Input:  vep_download_info -- [[:], vep_genome, vep_species, vep_cache_version]
		// Output: .cache -- downloaded and extracted VEP cache directory
		//

		if (download_vep_cache) {
			// Build input tuple: empty meta + genome/species/version for ENSEMBLVEP_DOWNLOAD
			vep_download_info = channel.of([[],vep_genome, vep_species, vep_cache_version])
			// vep_download_info: [[:], genome_str, species_str, cache_version_int]

			ENSEMBLVEP_DOWNLOAD (
				vep_download_info
			)

			ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
			ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

		}
		else {
			// Validate that the local cache directory exists and resolve the correct subdirectory
			// The annotation-cache S3 bucket uses a version-prefixed path; local paths do not
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

			// Collect the resolved cache root as a channel value
			ensemblvep_cache = channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
		}
		// ensemblvep_cache: path (or list-of-paths) to the VEP cache root directory

	emit:
		vep_cache = ensemblvep_cache  // path  -- VEP cache directory (downloaded or validated local)
		versions = ch_versions

}
