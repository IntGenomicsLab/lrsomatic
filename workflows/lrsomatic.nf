/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_lrsomatic_pipeline'
include { getGenomeAttribute     } from '../subworkflows/local/utils_nfcore_lrsomatic_pipeline'

//
// IMPORT MODULES
//
include { SAMTOOLS_CAT                      } from '../modules/nf-core/samtools/cat/main'
include { MINIMAP2_INDEX                    } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                    } from '../modules/nf-core/minimap2/align/main'
include { CRAMINO as CRAMINO_PRE            } from '../modules/local/cramino/main'
include { CRAMINO as CRAMINO_POST           } from '../modules/local/cramino/main'
include { NANOPLOT as NANOPLOT_PRE          } from '../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_POST         } from '../modules/nf-core/nanoplot/main'
include { MOSDEPTH                          } from '../modules/nf-core/mosdepth/main'
include { ASCAT                             } from '../modules/nf-core/ascat/main'
include { SEVERUS                           } from '../modules/nf-core/severus/main.nf'
include { METAEXTRACT                       } from '../modules/local/metaextract/main'
include { WAKHAN                            } from '../modules/local/wakhan/main'
include { FIBERTOOLSRS_PREDICTM6A           } from '../modules/local/fibertoolsrs/predictm6a'
include { FIBERTOOLSRS_FIRE                 } from '../modules/local/fibertoolsrs/fire'
include { FIBERTOOLSRS_NUCLEOSOMES          } from '../modules/local/fibertoolsrs/nucleosomes'
include { FIBERTOOLSRS_QC                   } from '../modules/local/fibertoolsrs/qc'
include { ENSEMBLVEP_VEP as SOMATIC_VEP     } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as GERMLINE_VEP    } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as SV_VEP          } from '../modules/nf-core/ensemblvep/vep/main.nf'
//
// IMPORT SUBWORKFLOWS
//
include { PREPARE_REFERENCE_FILES   } from '../subworkflows/local/prepare_reference_files'
include { PREPARE_ANNOTATION        } from '../subworkflows/local/prepare_annotation'
include { BAM_STATS_SAMTOOLS        } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { TUMOR_NORMAL_HAPPHASE     } from '../subworkflows/local/tumor_normal_happhase'
include { TUMOR_ONLY_HAPPHASE       } from '../subworkflows/local/tumor_only_happhase'





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LRSOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    // Channel format is [[meta], [bam]].
    // Where [meta] is [id, paired_data, method, specs, type]

    main:

    def clair3_modelMap = [
        'dna_r10.4.1_e8.2_400bps_sup@v5.2.0': 'r1041_e82_400bps_sup_v520',
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'r1041_e82_400bps_sup_v500',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'r1041_e82_400bps_sup_v430',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'r1041_e82_400bps_sup_v420',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'r1041_e82_400bps_sup_v410',
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'r1041_e82_260bps_sup_v400',
        'hifi_revio'                        : 'hifi'
    ]

    def clairs_modelMap = [
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v5.2.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'hifi_revio'                        : 'hifi_revio_ssrs'

    ]

    // Load in igenomes
    params.fasta = getGenomeAttribute('fasta')
    params.genome_name = getGenomeAttribute('genome_name')
    params.ascat_allele_files = getGenomeAttribute('ascat_alleles')
    params.ascat_loci_files = getGenomeAttribute('ascat_loci')
    params.ascat_gc_file = getGenomeAttribute('ascat_loci_gc')
    params.ascat_rt_file = getGenomeAttribute('ascat_loci_rt')
    params.centromere_bed = getGenomeAttribute('centromere_bed')
    params.pon_file = getGenomeAttribute('pon_file')
    params.bed_file = getGenomeAttribute('bed_file')
    params.vep_genome = getGenomeAttribute('vep_genome')
    params.vep_species = getGenomeAttribute('vep_species')
    params.dbsnp = getGenomeAttribute('dbsnp')
    params.colors = getGenomeAttribute('colors')
    params.onekgenomes = getGenomeAttribute('onekgenomes')
    params.gnomad = getGenomeAttribute('gnomad')

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: METAEXTRACT
    //
    // extracts the base calling model from the bam files

    METAEXTRACT( ch_samplesheet )

    basecall_meta = METAEXTRACT.out.meta_ext
    // [meta, basecall_model_str, kinetics_str]  -- basecall model and kinetics extracted from BAM header
    // Adds the base calling model to meta.basecall_model

    ch_samplesheet
        .join(basecall_meta)
        .map { meta, bam, basecall_model_meta, kinetics_meta ->
            def chosen_clair3_model = meta.clair3_model ?: clair3_modelMap.get(basecall_model_meta)
            def chosen_clairSTO_model = meta.clairSTO_model ?: clairs_modelMap.get(basecall_model_meta)
            def chosen_clairS_model = meta.clairS_model ?: clairs_modelMap.get(basecall_model_meta)
            def meta_new =[ id: meta.id,
                            paired_data: meta.paired_data,
                            type: meta.type,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            replicate: meta.replicate,
                            clair3_model: chosen_clair3_model,
                            clairS_model: chosen_clairS_model,
                            clairSTO_model: chosen_clairSTO_model,
                            kinetics: kinetics_meta]
            return[ meta_new, bam ]
        }
        .groupTuple()
        .map { meta, bam ->
            [ meta, bam.flatten()]
            }
        .set{ch_samplesheet}
    // [meta_full, [bam...]]  -- meta now includes: id, paired_data, type, platform, sex, fiber, clair3_model, clairS_model, clairSTO_model, kinetics



    //
    // SUBWORKFLOW: PREPARE_REFERENCE_FILES
    //

    PREPARE_REFERENCE_FILES (
        params.fasta,
        params.ascat_allele_files,
        params.ascat_loci_files,
        params.ascat_gc_file,
        params.ascat_rt_file,
        basecall_meta,
        clair3_modelMap
    )

    downloaded_clair3_models = PREPARE_REFERENCE_FILES.out.downloaded_clair3_models

    if (!params.skip_qc && !params.skip_cramino) {
        CRAMINO_PRE( ch_samplesheet )
        NANOPLOT_PRE(CRAMINO_PRE.out.arrow)
    }

    ch_samplesheet
        .map{ meta, bam ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'type',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, bam]
        }
        .set{ch_samplesheet_no_rep}


    // ch_samplesheet -> meta: [id, paired_data, platform, sex, type, fiber, basecall_model]
    //                   bam:  list of unaligned bams

    ch_split = ch_samplesheet_no_rep
        .branch { _meta, bam ->
            single: bam.size() == 1
            multiple: bam.size() > 1
        }

    //
    // MODULE: SAMTOOLS_CAT
    //
    // concatenates bam files from single sample

    SAMTOOLS_CAT ( ch_split.multiple )
        .bam
        .mix ( ch_split.single )
        .set { ch_cat_ubams }
    // [meta, bam]  -- single merged unaligned BAM per sample

    vep_cache = channel.empty()

    if (!params.skip_vep) {

        channel
            .of([
                vep_cache:          params.vep_cache,
                vep_cache_version:  params.vep_cache_version,
                vep_genome:         params.vep_genome,
                vep_args:           params.vep_args,
                vep_species:        params.vep_species,
                download_vep_cache: params.download_vep_cache
            ])

        PREPARE_ANNOTATION (
            params.vep_cache,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_args,
            params.vep_species,
            params.download_vep_cache
        )
        ch_versions = ch_versions.mix(PREPARE_ANNOTATION.out.versions)
        vep_cache = PREPARE_ANNOTATION.out.vep_cache.map {cache -> [[:], cache] }

    }

    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)
    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    ch_fai = PREPARE_REFERENCE_FILES.out.prepped_fai

    // ASCAT files
    allele_files = PREPARE_REFERENCE_FILES.out.allele_files
    loci_files = PREPARE_REFERENCE_FILES.out.loci_files
    gc_file = PREPARE_REFERENCE_FILES.out.gc_file
    rt_file = PREPARE_REFERENCE_FILES.out.rt_file

    //
    // MODULE: FIBERTOOLSRS_PREDICTM6A
    //
    // predict m6a in unaligned bam

    if (!params.skip_fiber) {
        if(!params.normal_fiber){
            ch_cat_ubams
            .branch { meta, _bams ->
                normal: meta.type == "normal"
                tumor: meta.type == "tumor"
                }
            .set { ch_cat_ubams_normal_branching }

            normal_bams = ch_cat_ubams_normal_branching.normal
            ubams = ch_cat_ubams_normal_branching.tumor
        }
        else {
            ubams = ch_cat_ubams
        }
            ubams
            .branch{ meta, _bams ->
                pacBio: meta.platform == "pb"
                ont: meta.platform == "ont"
            }
            .set{ch_cat_ubams_pacbio_ont_branching}

        pacbio_bams = ch_cat_ubams_pacbio_ont_branching.pacBio
        pacbio_bams
            .branch{meta, _bams ->
                kinetics: meta.kinetics == "true"
                noKinetics: meta.kinetics == "false"
            }
            .set{pacbio_bams}

        if (!params.skip_m6a) {
            FIBERTOOLSRS_PREDICTM6A (
                pacbio_bams.kinetics
            )
            pacbio_bams.noKinetics
                .mix(FIBERTOOLSRS_PREDICTM6A.out.bam)
                .set{predicted_bams}
        }
        else {
            pacbio_bams.noKinetics
                .mix(pacbio_bams.kinetics)
                .set{predicted_bams}
        }



        ch_cat_ubams_pacbio_ont_branching.ont
            .mix(predicted_bams)
            .set{fiber_branch}

        fiber_branch
            .branch{ meta, _bams ->
                fiber: meta.fiber == "y"
                nonFiber: meta.fiber == "n"
            }
            .set{fiber_branch}

        //
        // MODULE: FIBERTOOLSRS_NUCLEOSOMES
        //

        FIBERTOOLSRS_NUCLEOSOMES (
            fiber_branch.fiber
        )

        //
        // MODULE: FIBERTOOLSRS_FIRE
        //

        FIBERTOOLSRS_FIRE (
            FIBERTOOLSRS_NUCLEOSOMES.out.bam
        )

        if(!params.normal_fiber){
            fiber_branch.nonFiber
            .mix(normal_bams)
            .mix(FIBERTOOLSRS_FIRE.out.bam)
            .set{ch_cat_ubams}

        }
        else {
            fiber_branch.nonFiber
            .mix(FIBERTOOLSRS_FIRE.out.bam)
            .set{ch_cat_ubams}

        }

        if(!params.skip_qc) {
            //
            // MODULE: FIBERTOOLSRS_QC
            //

            FIBERTOOLSRS_QC (
                FIBERTOOLSRS_FIRE.out.bam
            )
        }
    }
    //
    // MODULE: MINIMAP2_ALIGN
    //
    // Aligns ubams
    // ch_cat_ubams: [meta, bam]  -- may include m6A/nucleosome/FIRE annotations for fiber-seq samples

    MINIMAP2_ALIGN (
        ch_cat_ubams,
        ch_fasta,
        true,
        'bai',
        "",
        ""
    )
    MINIMAP2_ALIGN.out.bam
        .set { ch_minimap_bam }
    // [meta, bam]  -- aligned BAM

    // ch_minimap_bams into tumor and paired to phase the paired ones on normal
    // and add index

    ch_minimap_bam
        .join(MINIMAP2_ALIGN.out.index)
        .branch { meta, _bams, _bais ->
                paired: meta.paired_data
                tumor_only: !meta.paired_data
        }
        .set { branched_minimap }
    // branched_minimap.paired:     [meta, bam, bai]  -- one item per sample (tumor AND normal flow separately)
    // branched_minimap.tumor_only: [meta, bam, bai]

    //
    // SUBWORFKLOW: TUMOR_NORMAL_HAPPHASE
    //
    // Phasing/haplotaging/small germline variant calling for tumor-normal samples

    TUMOR_NORMAL_HAPPHASE (
        branched_minimap.paired,
        ch_fasta,
        ch_fai,
        downloaded_clair3_models
    )

    ch_versions = ch_versions.mix(TUMOR_NORMAL_HAPPHASE.out.versions)

    //
    // SUBWORKFLOW: TUMOR_ONLY_HAPPHASE
    //
    // Phasing/haplotagging for tumor only samples

    dbsnp = file(params.dbsnp)
    colors = file(params.colors)
    onekgenomes = file(params.onekgenomes)
    gnomad = file(params.gnomad)

    TUMOR_ONLY_HAPPHASE (
        branched_minimap.tumor_only,
        ch_fasta,
        ch_fai,
        dbsnp,
        colors,
        onekgenomes,
        gnomad
    )

    germline_vep = TUMOR_NORMAL_HAPPHASE.out.germline_vep.mix(TUMOR_ONLY_HAPPHASE.out.germline_vep)
    // [meta, vcf, []]  -- germline variants merged from T/N and tumor-only paths
    somatic_vep = TUMOR_NORMAL_HAPPHASE.out.somatic_vep.mix(TUMOR_ONLY_HAPPHASE.out.somatic_vep)
    // [meta, vcf, []]  -- somatic variants merged from T/N and tumor-only paths

    if (!params.skip_vep) {
        //
        // MODULE: GERMLINE_VEP
        //
        if (params.vep_custom != null) {
            vep_custom = file(params.vep_custom)
        } else {
            vep_custom = []
        }
        if (params.vep_custom_tbi != null) {
            vep_custom_tbi = file(params.vep_custom_tbi)
        } else {
            vep_custom_tbi = []
        }
        GERMLINE_VEP (
            germline_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )

        //
        // MODULE: SOMATIC_VEP
        //

        SOMATIC_VEP (
            somatic_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )
    }

    ch_versions = ch_versions.mix(TUMOR_ONLY_HAPPHASE.out.versions)

    // Get Severus input channel
    TUMOR_NORMAL_HAPPHASE.out.tumor_normal_severus
        .mix(TUMOR_ONLY_HAPPHASE.out.tumor_only_severus)
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf, tbi ->
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf, tbi]
        }
        .set { severus_reformat }
    // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, phased_tbi]  -- normal_bam/bai are [] for tumor-only

    //
    // MODULE: SEVERUS
    //

    SEVERUS (
        severus_reformat,
        [[:], params.bed_file, params.pon_file]
    )

    ch_versions = ch_versions.mix(SEVERUS.out.versions)

    SEVERUS.out.all_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { sv_vep }
    // [meta, severus_all_vcf, []]  -- all SVs for VEP annotation

    if(!params.skip_vep) {
        SV_VEP (
            sv_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )
    }

    //
    // MODULE: CRAMINO
    //

    if (!params.skip_qc && !params.skip_cramino) {

        CRAMINO_POST ( ch_minimap_bam )
        NANOPLOT_POST(CRAMINO_POST.out.arrow)

    }

    //
    // Module: MOSDEPTH
    //

    ch_mosdepth_global = channel.empty()
    ch_mosdepth_summary = channel.empty()

    if (!params.skip_qc && !params.skip_mosdepth) {

        // prepare mosdepth input channel: we need to specify compulsory path to bed as well
        ch_minimap_bam.join(MINIMAP2_ALIGN.out.index)
            .map { meta, bam, bai -> [meta, bam, bai, []] }
            .set { ch_mosdepth_in }
        // [meta, bam, bai, []]  -- [] is the required empty BED path for MOSDEPTH

        MOSDEPTH (
            ch_mosdepth_in,
            ch_fasta
        )

        ch_mosdepth_global = MOSDEPTH.out.global_txt
        ch_mosdepth_summary = MOSDEPTH.out.summary_txt
    }

    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS
    //
    ch_bam_stats = channel.empty()
    ch_bam_flagstat = channel.empty()
    ch_bam_idxstats = channel.empty()

    if (!params.skip_qc && !params.skip_bamstats ) {

        BAM_STATS_SAMTOOLS (
            ch_minimap_bam.join(MINIMAP2_ALIGN.out.index), // Join bam channel with index channel
            ch_fasta
        )

        ch_bam_stats = BAM_STATS_SAMTOOLS.out.stats
        ch_bam_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
        ch_bam_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
    }

    //
    // MODULE: ASCAT
    //

    if (!params.skip_ascat) {
        severus_reformat
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _vcf ->
                return [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
            }
            .set { ascat_ch }
        // [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]  -- NOTE: normal before tumor (ASCAT convention)

        ASCAT (
            ascat_ch,
            params.genome_name,
            allele_files,
            loci_files,
            [],
            [],
            gc_file,
            rt_file
        )

        ch_versions = ch_versions.mix(ASCAT.out.versions)
    }

    //
    // MODULE: WAKHAN
    //

    if (!params.skip_wakhan) {

        // Prepare input channel for WAKHAN
        severus_reformat
            .join(SEVERUS.out.all_vcf)
            .set { wakhan_input }
        // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, phased_tbi, severus_all_vcf]

        WAKHAN (
            wakhan_input,
            ch_fasta,
            file(params.centromere_bed)
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'lrsomatic_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_stats.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_flagstat.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_idxstats.collect{it -> it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_global.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_summary.collect{it -> it[1]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files
            .collect()
            .map { files ->
                def multiqc_config_files = [file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)]
                if (params.multiqc_config) {
                    multiqc_config_files += [file(params.multiqc_config, checkIfExists: true)]
                }
                def multiqc_logo_file = params.multiqc_logo ? [file(params.multiqc_logo, checkIfExists: true)] : []
                [[id: 'multiqc'], files, multiqc_config_files, multiqc_logo_file, [], []]
            }
    )

    emit:
        multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html
        versions       = ch_versions                 // channel: [ path(versions.yml) ]



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
