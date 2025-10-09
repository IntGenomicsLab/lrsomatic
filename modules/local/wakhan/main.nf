// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process WAKHAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wakhan:0.2.0--pyhdfd78af_1':
        'biocontainers/wakhan:0.2.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(tumor_input), path(tumor_index), path(normal_input), path(normal_index), path(vcf), path(breakpoints)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*/*_genes_genome.html")                              , emit: genes_genome_html
    tuple val(meta), path("*/*_genes_genome.pdf")                               , emit: genes_genome_pdf
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints.html")            , emit: breakpoints_html
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints.pdf")             , emit: breakpoints_pdf
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints_subclonal.html")  , emit: breakpoints_subclonal_html
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints_subclonal.pdf")   , emit: breakpoints_subclonal_pdf
    tuple val(meta), path("*/*_genome_copynumbers_details.html")                , emit: copynumbers_details_html
    tuple val(meta), path("*/*_genome_copynumbers_details.pdf")                 , emit: copynumbers_details_pdf
    tuple val(meta), path("*/bed_output/*.bed")                                 , emit: bed_files
    tuple val(meta), path("*/variation_plots/*.html")                           , emit: variation_plots
    tuple val(meta), path("*/vcf_output/*_wakhan_cna_*.vcf")                    , emit: vcf_files
    tuple val(meta), path("*_heatmap_ploidy_purity.html")                       , emit: heatmap_html
    tuple val(meta), path("*_heatmap_ploidy_purity.html.pdf")                   , emit: heatmap_pdf
    tuple val(meta), path("*_optimized_peak.html")                              , emit: optimized_peak_html
    tuple val(meta), path("coverage_data/*.csv")                                , emit: coverage_csv
    tuple val(meta), path("coverage_plots/*.html")                              , emit: coverage_plots_html
    tuple val(meta), path("coverage_plots/*.pdf")                               , emit: coverage_plots_pdf
    tuple val(meta), path("phasing_output/*.html")                              , emit: phasing_html
    tuple val(meta), path("phasing_output/*.pdf")                               , emit: phasing_pdf
    tuple val(meta), path("phasing_output/*.rephased.vcf.gz")                   , emit: rephased_vcf
    tuple val(meta), path("phasing_output/*.rephased.vcf.gz.csi")               , emit: rephased_vcf_index
    tuple val(meta), path("snps_loh_plots/*_genome_snps_ratio_loh.html")        , emit: snps_loh_plot,      optional: true
    tuple val(meta), path("solutions_ranks.tsv")                                , emit: solutions_ranks
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def phased_vcf = normal_input ? "--normal-phased-vcf $vcf" : "--tumor-vcf $vcf"
    // WARN: Version information not provided by tool on CLI. Please update this string when upgrading BLAZE code
    def VERSION = "0.1.2"
    """
    wakhan \\
        --target-bam ${tumor_input} \\
        --breakpoints ${breakpoints} \\
        --reference ${reference} \\
        --genome-name ${prefix} \\
        --out-dir-plots . \\
        ${phased_vcf} \\
        ${args} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.1.2"

    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: $VERSION
    END_VERSIONS
    """
}
