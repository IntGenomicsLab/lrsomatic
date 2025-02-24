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
include { SAMTOOLS_CAT           } from '../modules/nf-core/samtools/cat/main'
include { MINIMAP2_INDEX         } from '../modules/nf-core/minimap2/index/main'
include { CLAIR3                 } from '../modules/local/clair3/main'
include { LONGPHASE_PHASE        } from '../modules/nf-core/longphase/phase/main'
include { LONGPHASE_HAPLOTAG     } from '../modules/nf-core/longphase/haplotag/main'
include { SEVERUS                } from '../modules/nf-core/longphase/haplotag/main'


include { MINIMAP2_ALIGN      } from '../modules/nf-core/minimap2/align/main'
include { CRAMINO as CRAMINO_PRE; CRAMINO as CRAMINO_POST } from '../modules/local/cramino/main'
include { MOSDEPTH         } from '../modules/nf-core/mosdepth/main'
include { METAEXTRACT         } from '../modules/local/metaextract/main'

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
   // METAEXTRACT(ch_samplesheet) | view { message -> "I say... $message" }
    basecall_meta = METAEXTRACT(ch_samplesheet)

    ch_samplesheet
    .join(basecall_meta)
    .map{ meta, bam, meta_ext ->
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
    */
    CRAMINO_PRE ( ch_cat_ubams )

    ch_versions = ch_versions.mix(CRAMINO_PRE.out.versions)
    
    
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

    
    // Need to join the tumor and normal information back together
    // since both bams need to be phased on the normal clair3 output
    // if it exists
    ch_minimap_bam
    .join(MINIMAP2_ALIGN.out.index) 
    .map{ meta, bam, bai ->
        def new_meta = [id: meta.id, 
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model
                        ]
        def bam_check = bam ? bam : []
        return[new_meta , [[type: meta.type], bam_check], [[type: meta.type], bai]]
    }
    .groupTuple()
    .branch {meta, bams, bais ->
        paired: meta.paired_data
        tumor_only : !meta.paired_data
    }
    .set{reformat_samples}

    // for the paired samples, only send the normal to CLAIR3
    reformat_samples.paired
    .map { meta, bams, bais ->
        def normal_bam = bams[0][0].type == "normal" ? bams[0][1] : bams[1][1]
        def normal_bai = bais[0][0].type == "normal" ? bais[0][1] : bais[1][1]
        return [meta, normal_bam, normal_bai]
    }
    .set{clair3_reformat_paired}

    // for the tumor only, send the tumor sample to CLAIR3
    reformat_samples.tumor_only
    .map{meta, bam, bai ->
        def tumor_bam = bam[0][1]
        def tumor_bai = bai[0][1]
        return [ meta, tumor_bam, tumor_bai]
    }
    .mix(clair3_reformat_paired)
    .set{clair3_reformat_samples}
    // FORMAT IS [ meta, normal_bam, normal_bai]
    // (except if its tumor only)

    CLAIR3 (
        clair3_reformat_samples,
        ch_fasta,
        ch_fai
    )

    // Now we need to join the germline snps to the tumor and normal samples
    reformat_samples.paired
    .map { meta, bams, bais ->
        def normal_bam = bams[0][0].type == "normal" ? bams[0][1] : bams[1][1]
        def tumor_bam = bams[0][0].type == "tumor" ? bams[0][1] : bams[1][1]
        def normal_bai = bais[0][0].type == "normal" ? bais[0][1] : bais[1][1]
        def tumor_bai = bais[0][0].type == "tumor" ? bais[0][1] : bais[1][1]

        return [ meta, normal_bam, tumor_bam, normal_bai, tumor_bai ]
    }
    .set{paired_samples}

    reformat_samples.tumor_only
    .map{meta, bam, bai ->
        def tumor_bam = bam[0][1]
        def tumor_bai = bai[0][1]
        return [ meta, [], tumor_bam, [], tumor_bai]
    }
    .mix(paired_samples)
    .join(CLAIR3.out.germline_vcf)
    .flatMap { meta, normal_bam, tumor_bam, normal_bai, tumor_bai, germline_vcf ->
            def meta_tumor = meta.clone()
            meta_tumor.type = 'tumor'
            def result = [[meta_tumor, tumor_bam, tumor_bai, germline_vcf]]
            
            if (normal_bam) {
                def meta_normal = meta.clone()
                meta_normal.type = 'normal'
                result << [meta_normal, normal_bam, normal_bai, germline_vcf]
            }
            return result
        } 
    .map{ meta, bam, bai, snps ->
        def snvs = []
        def mods = []
        return [ meta, bam, bai, snps, snvs, mods]
    }
    .set{longphase_reformat}
    // FORMAT IS [ meta, bam, bai, snps, [], [] ]

    LONGPHASE_PHASE(
        longphase_reformat,
        ch_fasta,
        ch_fai
    )

    longphase_reformat
    .join(LONGPHASE_PHASE.out.vcf)
    .map{meta, bam, bai, snps, snvs, mods, phased_snps ->
        return[meta, bam, bai,phased_snps, snvs, mods]
    }
    .set{longphase_reformat}
    // FORMAT IS [ meta, bam, bai, phased_snps, [], [] ]

    LONGPHASE_HAPLOTAG(
        longphase_reformat,
        ch_fasta,
        ch_fai
    )
    
    longphase_reformat
    .join(LONGPHASE_HAPLOTAG.out.bam)
    .map{meta, bam, bai, phased_snps, snvs, mods, hap_bam->
        return[meta, hap_bam, bai, phased_snps]
    }
    .map{ meta, hap_bam, bai, vcf ->
        def new_meta = [id: meta.id, 
                        paired_data: meta.paired_data,
                        platform: meta.platform,
                        basecall_model: meta.basecall_model
                        ]
        def bam_check = hap_bam ? hap_bam : []
        return[new_meta , [[type: meta.type], bam_check], [[type: meta.type], bai], vcf]
    }
    .groupTuple()
    .branch {meta, bams, bais, vcf->
        paired: meta.paired_data
        tumor_only : !meta.paired_data
    }
    .set{severus_reformat}

    severus_reformat.paired
    .map { meta, bams, bais, vcf ->
        def normal_bam = bams[0][0].type == "normal" ? bams[0][1] : bams[1][1]
        def tumor_bam = bams[0][0].type == "tumor" ? bams[0][1] : bams[1][1]
        def normal_bai = bais[0][0].type == "normal" ? bais[0][1] : bais[1][1]
        def tumor_bai = bais[0][0].type == "tumor" ? bais[0][1] : bais[1][1]

        return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf ]
    }
    .set{paired_samples}

    severus_reformat.tumor_only
    .map{meta, bam, bai, vcf->
        def tumor_bam = bam[0][1]
        def tumor_bai = bai[0][1]
        return [ meta,tumor_bam, tumor_bai, [], [], vcf]
    }
    .mix(paired_samples)
    .set{severus_reformat}
    // FORMAT IS [meta, tumor_hapbam, tumor_bai, normal_hapbam, normal_bai, vcf]



    
    // The channel is now [[meta], [bam]] With meta consisting of [id, paired_data, method, specs, type]
    
    // TODO: Add post-alignment QC step here, maybe add a subworkflow with all post-alignment QC together
    // 
    // MODULE: CRAMINO
    // 
    
    CRAMINO_POST ( ch_minimap_bam )

    ch_versions = ch_versions.mix(CRAMINO_POST.out.versions)


    
    
    //
    // Module: MOSDEPTH
    //
    
    // prepare mosdepth input channel: we need to specify compulsory path to bed as well
    ch_minimap_bam.join(MINIMAP2_ALIGN.out.index).map{ meta, bam, bai -> [meta, bam, bai, []]}.set{ch_mosdepth_in}

    MOSDEPTH ( ch_mosdepth_in,
        ch_fasta )  

    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    
    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS
    //
    BAM_STATS_SAMTOOLS (
        ch_minimap_bam.join(MINIMAP2_ALIGN.out.index), // Join bam channel with index channel
        ch_fasta
    )
    
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    
    
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