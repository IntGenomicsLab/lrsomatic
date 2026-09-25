//
// Stage the AmpliconArchitect data repository AmpliconClassifier reads at runtime
//

include { WGET as WGET_AA_DATA_REPO } from '../../modules/nf-core/wget/main'
include { UNTAR as UNTAR_AA_DATA_REPO } from '../../modules/nf-core/untar/main'

workflow PREPARE_AA_DATA_REPO {

    take:
    data_repo // path or null -- params.aa_data_repo
    repo_url  // URL or null  -- genome attribute aa_data_repo_url

    main:
    ch_versions = channel.empty()

    if (data_repo) {
        ch_data_repo = channel.value([ [ id: 'aa_data_repo' ], file(data_repo, checkIfExists: true) ])
    }
    else if (repo_url) {
        //
        // MODULES: WGET_AA_DATA_REPO -> UNTAR_AA_DATA_REPO (labels: process_single)
        // ~1.1 GB tarball; the plain build, not GRCh38_indexed, whose extra BWA index AC never reads
        //
        WGET_AA_DATA_REPO (
            channel.value([ [ id: 'aa_data_repo' ], repo_url ])
        )

        UNTAR_AA_DATA_REPO (
            WGET_AA_DATA_REPO.out.outfile
        )

        // .first(): UNTAR emits a queue channel, and every sample's classifier task
        // needs the same repo -- without this only the first sample would get it.
        ch_data_repo = UNTAR_AA_DATA_REPO.out.untar.first()
        ch_versions = ch_versions.mix(WGET_AA_DATA_REPO.out.versions, UNTAR_AA_DATA_REPO.out.versions)
    }
    else {
        error("AmpliconClassifier needs an AmpliconArchitect data repository, which is not published for ${params.genome}. Set --aa_data_repo <path>, or use --skip_ampliconclassifier.")
    }

    emit:
    data_repo = ch_data_repo   // [[id:'aa_data_repo'], dir]
    versions  = ch_versions
}
