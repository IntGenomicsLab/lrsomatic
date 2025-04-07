/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_lr_somatic_pipeline'

//
// IMPORT MODULES
//
include { SAMTOOLS_CAT              } from '../modules/nf-core/samtools/cat/main'
include { MINIMAP2_INDEX            } from '../modules/nf-core/minimap2/index/main'
include { CLAIRSTO                  } from '../modules/local/clairsto/main'
include { CLAIRS                    } from '../modules/local/clairs/main.nf'
include { MINIMAP2_ALIGN            } from '../modules/nf-core/minimap2/align/main'
include { CRAMINO as CRAMINO_PRE    } from '../modules/local/cramino/main'
include { CRAMINO as CRAMINO_POST   } from '../modules/local/cramino/main'
include { MOSDEPTH                  } from '../modules/nf-core/mosdepth/main'
include { ASCAT                     } from '../modules/nf-core/ascat/main'
include { CLAIR3                    } from '../modules/local/clair3/main'
include { LONGPHASE_PHASE           } from '../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG        } from '../modules/nf-core/longphase/haplotag/main'
include { SEVERUS                   } from '../modules/nf-core/severus/main.nf'
include { METAEXTRACT               } from '../modules/local/metaextract/main'
include { SAMTOOLS_INDEX            } from '../modules/nf-core/samtools/index/main.nf'
include { WAKHAN                    } from '../modules/local/wakhan/main'

//
// IMPORT SUBWORKFLOWS
//
include { PREPARE_REFERENCE_FILES   } from '../subworkflows/local/prepare_reference_files'
include { BAM_STATS_SAMTOOLS        } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { TUMOR_NORMAL_HAPPHASE     } from '../subworkflows/local/tumor_normal_happhase'
include { TUMOR_ONLY_HAPPHASE       } from '../subworkflows/local/tumor_only_happhase'



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
    // MODULE: METAEXTRACT
    //
    
    METAEXTRACT( ch_samplesheet )
    
    ch_versions  = ch_versions.mix(METAEXTRACT.out.versions)
    basecall_meta = METAEXTRACT.out.basecall_model

    ch_samplesheet
        .join(basecall_meta)
        .map { meta, bam, meta_ext ->
            def meta_new = meta + [ basecall_model: meta_ext]
            return[ meta_new, bam ]
        }
        .groupTuple()
        .map { meta, bam ->
            [ meta, bam.flatten()]
            }
        .set{ch_samplesheet}
    
    ch_split = ch_samplesheet
        .branch { meta, bam -> 
            single: bam.size() == 1
            multiple: bam.size() > 1
        }
    
    //
    // MODULE: SAMTOOLS_CAT
    //
    
    SAMTOOLS_CAT ( ch_split.multiple )
        .bam
        .mix ( ch_split.single )
        .set { ch_cat_ubams }
        
    ch_versions = ch_versions.mix(SAMTOOLS_CAT.out.versions)
    
    //
    // MODULE: CRAMINO
    //
    
    if (!params.skip_qc && !params.skip_cramino) {
        
        CRAMINO_PRE ( ch_cat_ubams )

        ch_versions = ch_versions.mix(CRAMINO_PRE.out.versions)
    }

    
    //
    // SUBWORKFLOW: PREPARE_REFERENCE_FILES
    //
    
    PREPARE_REFERENCE_FILES ( 
        params.fasta,
        params.ascat_allele_files,
        params.ascat_loci_files,
        params.ascat_gc_file,
        params.ascat_rt_file 
    )
    
    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)
    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    ch_fai = PREPARE_REFERENCE_FILES.out.prepped_fai
    
    // ASCAT files   
    allele_files = PREPARE_REFERENCE_FILES.out.allele_files
    loci_files = PREPARE_REFERENCE_FILES.out.loci_files
    gc_file = PREPARE_REFERENCE_FILES.out.gc_file
    rt_file = PREPARE_REFERENCE_FILES.out.rt_file
    
    
    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)

    //
    // MODULE: MINIMAP2_ALIGN
    //
    
    MINIMAP2_ALIGN ( 
        ch_cat_ubams,
        ch_fasta,
        true,
        'bai',
        "", 
        "" 
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.bam 
        .set { ch_minimap_bam } 

    // Need to join the tumor and normal information back together
    // since both bams need to be phased on the normal clair3 output
    // if it exists
    ch_minimap_bam
        .join(MINIMAP2_ALIGN.out.index)
        .branch { meta, bams, bais ->
                paired: meta.paired_data
                tumor_only: !meta.paired_data
        }
        .set { branched_minimap }

    //
    // SUBWORFKLOW: TUMOR_NORMAL_HAPPHASE
    //
    
    TUMOR_NORMAL_HAPPHASE (
        branched_minimap.paired,
        ch_fasta,
        ch_fai
    )
    
    ch_versions = ch_versions.mix(TUMOR_NORMAL_HAPPHASE.out.versions)
    
    //
    // SUBWORKFLOW: TUMOR_ONLY_HAPPHASE
    //
    
    TUMOR_ONLY_HAPPHASE (
        branched_minimap.tumor_only,
        ch_fasta,
        ch_fai
    )
    
    ch_versions = ch_versions.mix(TUMOR_NORMAL_HAPPHASE.out.versions)
    
    // Get Severus input channel
    TUMOR_NORMAL_HAPPHASE.out.tumor_normal_severus
        .mix(TUMOR_ONLY_HAPPHASE.out.tumor_only_severus)
        .set { severus_reformat }
    // Format is [meta, tumor_hapbam, tumor_bai, normal_hapbam, normal_bai, vcf]
    
    // Get ClairS input channel
    TUMOR_NORMAL_HAPPHASE.out.tumor_normal_severus
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
            return[meta , normal_bam, normal_bai, tumor_bam, tumor_bai]
        }
        .set { clairs_input }
        
    //
    // MODULE: CLAIRS
    //    
    
    CLAIRS (
        clairs_input,
        ch_fasta,
        ch_fai
    )
    
    ch_versions = ch_versions.mix(CLAIRS.out.versions)
    
    //
    // MODULE: SEVERUS
    //
    
    if (!params.skip_severus) {
        
        SEVERUS (
            severus_reformat,
            [[:], params.bed_file, params.pon_file]
        )
        
        ch_versions = ch_versions.mix(SEVERUS.out.versions)
    }
    
    // 
    // MODULE: CRAMINO
    // 
    
    if (!params.skip_qc && !params.skip_cramino) {
        
        CRAMINO_POST ( ch_minimap_bam )

        ch_versions = ch_versions.mix(CRAMINO_POST.out.versions)
    }
    
    //
    // Module: MOSDEPTH
    //
    
    if (!params.skip_qc && params.skip_mosdepth) {
        
        // prepare mosdepth input channel: we need to specify compulsory path to bed as well
        ch_minimap_bam.join(MINIMAP2_ALIGN.out.index)
            .map { meta, bam, bai -> [meta, bam, bai, []] }
            .set { ch_mosdepth_in }

        MOSDEPTH ( 
            ch_mosdepth_in,
            ch_fasta 
        )  

        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    }

    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS
    //
    
    if (!params.skip_qc && params.skip_bamstats ) {
        
        BAM_STATS_SAMTOOLS (
            ch_minimap_bam.join(MINIMAP2_ALIGN.out.index), // Join bam channel with index channel
            ch_fasta
        )
        
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }
    
    //
    // MODULE: ASCAT
    //
    
    if (!skip_ascat) {
        severus_reformat
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ->
                return[meta , normal_bam, normal_bai, tumor_bam, tumor_bai]
            }
            .set { ascat_ch }
        
        ASCAT (
            ascat_ch,
            allele_files,
            loci_files,
            [],
            [],
            [],
            []
        )
        
        ch_versions = ch_versions.mix(ASCAT.out.versions)
    }
    
    /*
    //
    // MODULE: WAKHAN
    //
    
    if (!skip_wakhan) {
        
        WAKHAN (
            severus_reformat,
            ch_fasta
        )
        
        ch_versions = ch_versions.mix(WAKHAN.out.versions)
    }

    
    */
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
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]}.ifEmpty([]))
    
    
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