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

include { SAMTOOLS_CAT as SAMTOOLS_CAT_TUMOUR  } from '../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_CAT as SAMTOOLS_CAT_NORMAL  } from '../modules/nf-core/samtools/cat/main'
include { MINIMAP2_INDEX                       } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                       } from '../modules/nf-core/minimap2/align/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LR_SOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_samplesheet
        .branch{
            meta, tumour ->
                single: tumour.size() == 1
                    return [ meta, tumour.flatten(), normal.flatten() ] // TODO: this probably doesnt work
                multiple: tumour.size() > 1
                    return [ meta, tumour.flatten(), normal.flatten() ]
        }
        .set { ch_ubams }
    
    // TODO: Split workflow here in paired/tumour-only?
    
    //
    // MODULE: Combine bam files from the same sample (TUMOUR ubams)
    //
    // TODO: Ensure it only takes tumour bam here
    SAMTOOLS_CAT_TUMOUR ( ch_ubams.multiple )
        .reads
        .mix ( ch_ubams.single )
        .set { ch_cat_ubams }

    ch_versions = ch_versions.mix (SAMTOOLS_CAT_TUMOUR.out.versions.first().ifEmpty(null))
    
    //
    // MODULE: Combine bam files from the same sample (NORMAL ubams)
    //
    // TODO: Ensure it only takes normal bam here
    SAMTOOLS_CAT_NORMAL ( ch_ubams.multiple )
        .reads
        .mix ( ch_ubams.single )
        .set { ch_cat_ubams }


    ch_versions = ch_versions.mix (SAMTOOLS_CAT_NORMAL.out.versions.first().ifEmpty(null))
    
    // TODO: Add pre-alignment QC step here
    //
    // MODULE: CRAMINO
    //
    CRAMINO_PRE ( )

    
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
    
    // Create minimap2 index channel
    
    
    if (!params.skip_save_minimap2_index) {
        
        MINIMAP2_INDEX ( ch_fasta )
        ch_minimap_index = MINIMAP2_INDEX.out.index
        
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MODULE: Run MINIMAP2_ALIGN
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
        | set { ch_minimap_bam }

    
    // TODO: Add post-alignment QC step here
    // 
    // MODULE: CRAMINO
    // 
    
    CRAMINO_POST ( )
    
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

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
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
