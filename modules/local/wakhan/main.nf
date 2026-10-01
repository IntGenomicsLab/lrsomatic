process WAKHAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // No biocontainer was built for 0.5.0; this image installs bioconda::wakhan=0.5.0
    // (recipe: https://github.com/ljwharbers/Wakhan/tree/container/0.5.0/container)
    container "ghcr.io/ljwharbers/wakhan:0.5.0"

    input:
    tuple val(meta), path(tumor_input), path(tumor_index), path(normal_input), path(normal_index), path(vcf), path(breakpoints)
    tuple val(meta2), path(reference)
    path(centromere_bed)

    output:
    // Per-solution outputs live in solution_<ploidy>_<purity>_<confidence>/; solution_rank_<n> links to them
    tuple val(meta), path("solution_*_*_*/integer_profile.html")                , emit: integer_profile_html
    tuple val(meta), path("solution_*_*_*/integer_profile.pdf")                 , emit: integer_profile_pdf
    tuple val(meta), path("solution_*_*_*/integer_profile.{bed,vcf}")           , emit: integer_profile
    tuple val(meta), path("solution_*_*_*/subclonal_profile.html")              , emit: subclonal_profile_html
    tuple val(meta), path("solution_*_*_*/subclonal_profile.pdf")               , emit: subclonal_profile_pdf
    tuple val(meta), path("solution_*_*_*/subclonal_profile.{bed,vcf}")         , emit: subclonal_profile
    tuple val(meta), path("solution_*_*_*/genes/genes_copynumber_states.*")     , emit: genes
    tuple val(meta), path("solution_*_*_*/HiScanner_plots_data.zip")            , emit: hiscanner_zip
    tuple val(meta), path("*_heatmap_ploidy_purity.html")                       , emit: heatmap_html
    tuple val(meta), path("*_heatmap_ploidy_purity.html.pdf")                   , emit: heatmap_pdf
    tuple val(meta), path("coverage_data/*.csv")                                , emit: coverage_csv
    tuple val(meta), path("coverage_data/*.png")                                , emit: coverage_png
    tuple val(meta), path("coverage_plots/*.html")                              , emit: coverage_plots_html
    tuple val(meta), path("coverage_plots/*.pdf")                               , emit: coverage_plots_pdf
    tuple val(meta), path("phasing_output/*.html")                              , emit: phasing_html
    tuple val(meta), path("phasing_output/*.pdf")                               , emit: phasing_pdf
    tuple val(meta), path("phasing_output/*rephased.vcf.gz")                    , emit: rephased_vcf
    tuple val(meta), path("phasing_output/*rephased.vcf.gz.csi")                , emit: rephased_vcf_index
    tuple val(meta), path("snps_loh_plots/*_genome_snps_ratio_loh.html")        , emit: snps_loh_plot,      optional: true
    tuple val(meta), path("solutions_ranks.tsv")                                , emit: solutions_ranks
    // Whole directories, not the plots inside: every solution's plot has the same basename,
    // and LRSOMATICREPORT resolves them by rank directory
    tuple val(meta), path("solution_rank_*", type: 'dir')                       , emit: solution_dirs,   optional: true
    tuple val("${task.process}"), val('wakhan'), eval("wakhan --version 2>/dev/null | sed 's/^wakhan //'"), topic: versions, emit: versions_wakhan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def phased_vcf = normal_input ? "--normal-phased-vcf $vcf" : "--tumor-phased-vcf $vcf"
    def centromere = centromere_bed ? "--centromere-bed \$PWD/${centromere_bed}" : ""

    """
    wakhan \\
        all \\
        --target-bam ${tumor_input} \\
        --breakpoints ${breakpoints} \\
        --reference ${reference} \\
        --genome-name ${prefix} \\
        --out-dir-plots . \\
        ${phased_vcf} \\
        ${centromere} \\
        ${args} \\
        --threads ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sol=solution_2.0_0.8_0.9
    mkdir -p \$sol/genes coverage_data coverage_plots phasing_output snps_loh_plots
    for f in integer_profile subclonal_profile; do
        touch \$sol/\$f.html \$sol/\$f.pdf \$sol/\$f.bed \$sol/\$f.vcf
    done
    touch \$sol/genes/genes_copynumber_states.html \$sol/genes/genes_copynumber_states.pdf \$sol/genes/genes_copynumber_states.bed
    touch \$sol/HiScanner_plots_data.zip
    ln -s \$sol solution_rank_1
    touch ${prefix}_heatmap_ploidy_purity.html ${prefix}_heatmap_ploidy_purity.html.pdf
    touch coverage_data/coverage.csv coverage_data/cn_peaks.png coverage_plots/COVERAGE_INDEX.html coverage_plots/chr1.pdf
    touch phasing_output/PHASE_CORRECTION_INDEX.html phasing_output/chr1.pdf
    echo "" | gzip > phasing_output/rephased.vcf.gz
    touch phasing_output/rephased.vcf.gz.csi
    touch snps_loh_plots/${prefix}_genome_snps_ratio_loh.html
    printf 'repository_name\\tdna_purity\\tcell_purity\\tploidy\\tconfidence\\tsolution_rank\\n%s\\t0.8\\t0.8\\t2.0\\t0.9\\t1\\n' \$sol > solutions_ranks.tsv
    """
}
