//
// Reshape the VEP plugin releases that ship as zip archives (REVEL and EVE), and download a remote
// ClinVar once so the VEP tasks never stage it from its host themselves
//

include { UNZIP as UNZIP_REVEL } from '../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_EVE   } from '../../modules/nf-core/unzip/main.nf'
include { WGET as WGET_REVEL   } from '../../modules/nf-core/wget/main'
include { WGET as WGET_EVE     } from '../../modules/nf-core/wget/main'
include { VEPPLUGIN_REVEL      } from '../../modules/local/vepplugin/revel/main.nf'
include { VEPPLUGIN_EVE        } from '../../modules/local/vepplugin/eve/main.nf'
include { VEPPLUGIN_CLINVAR    } from '../../modules/local/vepplugin/clinvar/main.nf'

workflow PREPARE_VEP_PLUGINS {

    take:
    plugins // map: the resolved plugin configuration from resolveVepPlugins()

    main:

    def prepare = plugins.prepare

    ch_versions = channel.empty()
    def staged = [ channel.fromList(plugins.ready_files) ]

    //
    // MODULES: WGET_REVEL -> UNZIP_REVEL -> VEPPLUGIN_REVEL (labels: process_single, process_single, process_medium)
    // Input:  the REVEL release, as a URL or as a local zip
    // Output: .files -- revel_grch38.tsv.gz and its index
    // Remote releases go through wget: REVEL 403s a request without a User-Agent and EVE redirects HTTPS to HTTP
    //
    if (prepare.containsKey('vep_revel')) {
        if (prepare['vep_revel'].toString().contains('://')) {
            WGET_REVEL (
                channel.value([ [ id: 'revel' ], prepare['vep_revel'] ])
            )

            ch_revel_zip = WGET_REVEL.out.outfile
            ch_versions = ch_versions.mix(WGET_REVEL.out.versions)
        }
        else {
            ch_revel_zip = channel.value([ [ id: 'revel' ], file(prepare['vep_revel'], checkIfExists: true) ])
        }

        UNZIP_REVEL (
            ch_revel_zip
        )

        VEPPLUGIN_REVEL (
            UNZIP_REVEL.out.unzipped_archive.map { _meta, dir -> dir }
        )

        staged << VEPPLUGIN_REVEL.out.files
        ch_versions = ch_versions.mix(UNZIP_REVEL.out.versions)
    }

    //
    // MODULES: WGET_EVE -> UNZIP_EVE -> VEPPLUGIN_EVE (labels: process_single, process_single, process_medium)
    // Input:  the EVE release, as a URL or as a local zip -- one VCF per protein
    // Output: .files -- eve_merged.vcf.gz and its index
    //
    if (prepare.containsKey('vep_eve')) {
        if (prepare['vep_eve'].toString().contains('://')) {
            WGET_EVE (
                channel.value([ [ id: 'eve' ], prepare['vep_eve'] ])
            )

            ch_eve_zip = WGET_EVE.out.outfile
            ch_versions = ch_versions.mix(WGET_EVE.out.versions)
        }
        else {
            ch_eve_zip = channel.value([ [ id: 'eve' ], file(prepare['vep_eve'], checkIfExists: true) ])
        }

        UNZIP_EVE (
            ch_eve_zip
        )

        VEPPLUGIN_EVE (
            UNZIP_EVE.out.unzipped_archive.map { _meta, dir -> dir }
        )

        staged << VEPPLUGIN_EVE.out.files
        ch_versions = ch_versions.mix(UNZIP_EVE.out.versions)
    }

    //
    // MODULE: VEPPLUGIN_CLINVAR (label: process_single)
    // Input:  the ClinVar VCF and index URLs, and the MD5s pinning them (either may be null)
    // Output: .files -- the VCF and its index, under the VCF's own basename
    // One download per run: a foreign file is re-checked on its host by GERMLINE_VEP and SOMATIC_VEP
    // for every sample, and NCBI answers the burst a multi-sample run sends with 503
    //
    if (prepare.containsKey('vep_clinvar')) {
        def clinvar = prepare['vep_clinvar']

        VEPPLUGIN_CLINVAR (
            channel.value([ clinvar.vcf, clinvar.tbi, clinvar.md5, clinvar.tbi_md5 ])
        )

        staged << VEPPLUGIN_CLINVAR.out.files
        ch_versions = ch_versions.mix(VEPPLUGIN_CLINVAR.out.versions)
    }

    // Value channel read by both VEP tasks; ifEmpty carries the no-plugins case, since collect() emits nothing then
    ch_extra_files = staged
        .inject(channel.empty()) { acc, ch -> acc.mix(ch) }
        .flatten()
        .collect()
        .ifEmpty([])

    emit:
    extra_files = ch_extra_files // channel: value list of plugin .pm and data files
    versions    = ch_versions    // channel: versions.yml files
}
