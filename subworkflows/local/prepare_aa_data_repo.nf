//
// Stage the AmpliconArchitect data repository AmpliconClassifier reads at runtime
//

include { AADATAREPO_DOWNLOAD         } from '../../modules/local/aadatarepo/download/main'
include { UNTAR as UNTAR_AA_DATA_REPO } from '../../modules/nf-core/untar/main'

workflow PREPARE_AA_DATA_REPO {

    take:
    data_repo // path or null -- params.aa_data_repo
    repo_url  // URL or null  -- genome attribute aa_data_repo_url
    repo_md5  // MD5 or null  -- --aa_data_repo_md5, else the genome attribute

    main:
    ch_versions = channel.empty()

    // The MD5 is checked by the download task, so a local repo would silently skip it
    if (data_repo && params.aa_data_repo_md5) {
        error("--aa_data_repo_md5: only checks a repository the pipeline downloads. Drop it when --aa_data_repo is set.")
    }

    if (data_repo) {
        ch_data_repo = channel.value([ [ id: 'aa_data_repo' ], file(data_repo, checkIfExists: true) ])
    }
    else if (repo_url) {
        //
        // MODULES: AADATAREPO_DOWNLOAD -> UNTAR_AA_DATA_REPO (labels: process_single)
        // ~1.1 GB tarball; the plain build, not GRCh38_indexed, whose extra BWA index AC never reads
        //
        if (!repo_md5) {
            log.warn("AmpliconClassifier: the data repository '${repo_url}' is downloaded without --aa_data_repo_md5, so the repository is not verified.")
        }
        AADATAREPO_DOWNLOAD (
            channel.value([ [ id: 'aa_data_repo' ], repo_url, repo_md5 ])
        )

        UNTAR_AA_DATA_REPO (
            AADATAREPO_DOWNLOAD.out.archive
        )

        // .first(): UNTAR emits a queue channel, and every sample's classifier task
        // needs the same repo -- without this only the first sample would get it.
        ch_data_repo = UNTAR_AA_DATA_REPO.out.untar.first()
        // UNTAR reports its version through the versions topic
        ch_versions = ch_versions.mix(AADATAREPO_DOWNLOAD.out.versions)
    }
    else {
        // An empty repo runs no classifier tasks; CoRAL reconstruction is unaffected
        log.warn("No AmpliconArchitect data repository is published for ${params.genome}: skipping AmpliconClassifier. Set --aa_data_repo <path> to classify.")
        ch_data_repo = channel.empty()
    }

    emit:
    data_repo = ch_data_repo   // [[id:'aa_data_repo'], dir], or empty when no repo is available
    versions  = ch_versions
}
