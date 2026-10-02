process SEVERUS_WHITELIST {
    tag "$meta.id"
    label 'process_medium'

    // No conda: the --whitelist fixes exist only in the patched image, and bioconda's severus would silently drop them
    // Severus 1.7 + github.com/AmberVerhasselt/Severus/tree/whitelist-reciprocal-corroboration (6813dee); revert to the biocontainer once KolmogorovLab/Severus carries it
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/amberverhasselt/severus-sif:1.7-whitelist-6813dee'
        : 'docker.io/amberverhasselt/severus:1.7-whitelist-6813dee'}"

    input:
    tuple val(meta), path(target_input), path(target_index), path(control_input), path(control_index), path(vcf), path(tbi)
    tuple val(meta2), path(bed), path(pon_path)
    tuple val(meta3), path(whitelist)

    output:
    tuple val(meta), path("${prefix}/severus.log")                              , emit: log
    tuple val(meta), path("${prefix}/read_qual.txt")                            , emit: read_qual
    tuple val(meta), path("${prefix}/breakpoints_double.csv")                   , emit: breakpoints_double
    tuple val(meta), path("${prefix}/read_alignments")                          , emit: read_alignments                  , optional: true
    tuple val(meta), path("${prefix}/read_ids.csv")                             , emit: read_ids                         , optional: true
    tuple val(meta), path("${prefix}/severus_collaped_dup.bed")                 , emit: collapsed_dup                    , optional: true
    tuple val(meta), path("${prefix}/severus_LOH.bed")                          , emit: loh                              , optional: true
    tuple val(meta), path("${prefix}/all_SVs/severus_all.vcf.gz")               , emit: all_vcf                          , optional: true
    tuple val(meta), path("${prefix}/all_SVs/breakpoint_clusters_list.tsv")    , emit: all_breakpoints_clusters_list    , optional: true
    tuple val(meta), path("${prefix}/all_SVs/breakpoint_clusters.tsv")         , emit: all_breakpoints_clusters         , optional: true
    tuple val(meta), path("${prefix}/all_SVs/plots/severus*.html")              , emit: all_plots                        , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/severus_somatic.vcf.gz")       , emit: somatic_vcf                      , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/breakpoint_clusters_list.tsv"), emit: somatic_breakpoints_clusters_list, optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/breakpoint_clusters.tsv")     , emit: somatic_breakpoints_clusters     , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/plots/severus*.html")          , emit: somatic_plots                    , optional: true
    // The -whitelist suffix keeps the patched build distinguishable from stock severus in the versions report
    tuple val("${task.process}"), val('severus'), eval('echo "$(severus --version)-whitelist-6813dee"'), emit: versions_severus, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SEVERUS_WHITELIST does not support Conda. Please use Docker / Singularity / Apptainer instead, or run without --severus_whitelist."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def control = control_input ? "--control-bam ${control_input}" : ""
    def vntr_bed = bed ? "--vntr-bed ${bed}" : ""
    def phasing_vcf = vcf ? "--phasing-vcf ${vcf}" : ""
    def pon = pon_path && (!control_input) ? "--PON ${pon_path}" : ""

    """
    severus \\
        $args \\
        --threads $task.cpus \\
        --target-bam $target_input \\
        $vntr_bed \\
        $pon \\
        $control \\
        $phasing_vcf \\
        --whitelist $whitelist \\
        --out-dir ${prefix}

    bgzip ${prefix}/somatic_SVs/severus_somatic.vcf
    tabix -p vcf ${prefix}/somatic_SVs/severus_somatic.vcf.gz
    bgzip ${prefix}/all_SVs/severus_all.vcf
    tabix -p vcf ${prefix}/all_SVs/severus_all.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/all_SVs/plots
    mkdir -p ${prefix}/somatic_SVs/plots

    touch ${prefix}/severus_collaped_dup.bed
    touch ${prefix}/severus.log
    touch ${prefix}/severus_LOH.bed
    touch ${prefix}/read_ids.csv
    touch ${prefix}/read_qual.txt
    touch ${prefix}/breakpoints_double.csv
    echo "" | gzip -n > ${prefix}/all_SVs/severus_all.vcf.gz
    touch ${prefix}/all_SVs/breakpoint_clusters_list.tsv
    touch ${prefix}/all_SVs/breakpoint_clusters.tsv
    touch ${prefix}/all_SVs/plots/severus_0.html
    echo "" | gzip -n > ${prefix}/somatic_SVs/severus_somatic.vcf.gz
    touch ${prefix}/somatic_SVs/breakpoint_clusters_list.tsv
    touch ${prefix}/somatic_SVs/breakpoint_clusters.tsv
    touch ${prefix}/somatic_SVs/plots/severus_0.html
    """
}
