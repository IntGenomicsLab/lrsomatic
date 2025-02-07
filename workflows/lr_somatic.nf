/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_lr_somatic_pipeline'

//
// IMPORT MODULES
//
include { SAMTOOLS_CAT        } from '../modules/nf-core/samtools/cat/main'
include { MINIMAP2_INDEX      } from '../modules/nf-core/minimap2/index/main'
include { CLAIRSTO            } from '../modules/local/clairsto/main'
include { MINIMAP2_ALIGN      } from '../modules/nf-core/minimap2/align/main'
include { CRAMINO as CRAMINO_PRE; CRAMINO as CRAMINO_POST       } from '../modules/local/cramino/main'
//
// IMPORT SUBWORKFLOWS
//
include { PREPARE_REFERENCE_FILES     } from '../subworkflows/local/prepare_reference_files'
include { BAM_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LR_SOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    // Channel format is [[meta], [bam]]. 
    // Where [meta] is [id, paired_data, method, specs, type]
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    

    //
    // MODULE: Combine bam files from the same sample 
    //
    
    // Take channels where there are multiple bam files in the list
    
    ch_split = ch_samplesheet
        .branch { meta, bam -> 
            single: bam.size() == 1
            multiple: bam.size() > 1
        }
    
    SAMTOOLS_CAT ( ch_split.multiple )
        .bam
        .mix ( ch_split.single )
        .set { ch_cat_ubams }
    ch_versions = ch_versions.mix (SAMTOOLS_CAT.out.versions)
    
    /*
    // TODO: Add pre-alignment QC step here
    // Maybe add a subworkflow with all pre-alignment QC together if there will be more than CRAMINO
    //
    // MODULE: CRAMINO
    //
    */
    CRAMINO_PRE ( ch_cat_ubams )
    
    
    //
    // SUBWORKFLOW: PREPARE_REFERENCE_FILES
    //
    
    PREPARE_REFERENCE_FILES ( params.fasta )
    
    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    ch_fai = PREPARE_REFERENCE_FILES.out.prepped_fai
    
    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)
    
    //
    // MODULE: Run MINIMAP2_INDEX
    //
    if (!params.skip_save_minimap2_index) {
        
        MINIMAP2_INDEX ( ch_fasta )
        ch_minimap_index = MINIMAP2_INDEX.out.index
        
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MODULE: MINIMAP2_ALIGN
    //
    MINIMAP2_ALIGN ( 
        ch_cat_ubams,
        ch_minimap_index,
        true,
        'bai',
        "", 
        "" 
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.bam 
        .set { ch_minimap_bam } 


    CLAIRSTO (
        ch_minimap_bam.join(MINIMAP2_ALIGN.out.index),
        ch_fasta,
        ch_fai
    )
    // The channel is now [[meta], [bam]] With meta consisting of [id, paired_data, method, specs, type]
    
    // TODO: Add post-alignment QC step here, maybe add a subworkflow with all post-alignment QC together
    // 
    // MODULE: CRAMINO
    // 
    
    //CRAMINO_POST ( )
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'lr_somatic_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files    = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    
    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]


    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/