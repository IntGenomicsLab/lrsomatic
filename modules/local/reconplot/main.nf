process RECONPLOT {
    tag "${meta.id}:${cn_source}+${sv_source}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Built from containers/reconplot/Dockerfile: R deps + ReConPlot package (not on conda; the wrapper is staged as source).
    // Override per site with `process { withName: '.*:RECONPLOT_(ASCAT_SEVERUS|WAKHAN_SEVERUS|SAVANA)' { container = ... } }`.
    container "ghcr.io/tim-yu/reconplot@sha256:1145fc5aebe0227bec371f4c59b08b9a09871498e403c01b83f83973149ae9e7"

    input:
    // cn_files / sv_files are the caller output files the wrapper's parsers discover by name.
    // When cn_source == sv_source (e.g. savana) put everything in cn_files and leave sv_files empty.
    tuple val(meta), val(cn_source), path(cn_files, stageAs: 'cn_input/*'), val(sv_source), path(sv_files, stageAs: 'sv_input/*')
    tuple val(meta2), path(reconplot_src)   // ReConPlot wrapper (contains run_reconplot.R + R/)
    tuple val(meta3), path(reconplot_pkg)   // ReConPlot R package source; installed only if the env lacks it
    val(genome)                             // hg38 | hg19 | T2T | mm10 | mm39

    output:
    tuple val(meta), path("${prefix}/per_chromosome/*.{pdf,png}"), emit: per_chromosome
    tuple val(meta), path("${prefix}/genome_wide/*.{pdf,png}")   , emit: genome_wide
    tuple val(meta), path("${prefix}/focus/*.{pdf,png}")         , emit: focus  , optional: true
    tuple val(meta), path("${prefix}/*.reconplot_{cn,sv}.tsv")   , emit: tables
    tuple val(meta), path("${prefix}/reconplot.log")             , emit: log
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''   // shared filters (e.g. --min-svlen, --max-cn)
    def args2  = task.ext.args2 ?: ''   // genome-wide strip extras
    def args3  = task.ext.args3 ?: ''   // focus panel extras (--regions/--genes/--baf-track); focus skipped if empty
    prefix     = task.ext.prefix ?: (cn_source == sv_source ? "${cn_source}" : "${cn_source}_${sv_source}")
    def sample = meta.id
    def source_args = cn_source == sv_source
        ? "--source ${cn_source} --input cn_input"
        : "--cn-source ${cn_source} --cn-input cn_input --sv-source ${sv_source} --sv-input sv_input"
    // Wakhan's parser expects <dir>/solutions_ranks.tsv + <dir>/<solution>/bed_output/*.bed
    def layout_cmd = cn_source == 'wakhan'
        ? "mkdir -p cn_input/solution_1/bed_output && mv cn_input/*.bed cn_input/solution_1/bed_output/"
        : ""
    def focus_cmd = args3
        ? """
    Rscript ${reconplot_src}/run_reconplot.R ${source_args} --sample ${sample} --prefix ${sample} \\
        --genome ${genome} --outdir ${prefix}/focus --layout together ${args} ${args3} 2>&1 | tee -a ${prefix}/reconplot.log
    """
        : ""

    """
    # ReConPlot is pre-installed in the container; conda envs get it from the staged source tree
    if ! Rscript -e 'suppressMessages(library(ReConPlot))' 2>/dev/null; then
        mkdir -p rlib
        R CMD INSTALL --no-docs --no-html -l rlib ${reconplot_pkg} > rlib_install.log 2>&1
        export R_LIBS=\$PWD/rlib\${R_LIBS:+:\$R_LIBS}
    fi

    ${layout_cmd}
    mkdir -p ${prefix}

    Rscript ${reconplot_src}/run_reconplot.R ${source_args} --sample ${sample} --prefix ${sample} \\
        --genome ${genome} --outdir ${prefix}/per_chromosome --regions all --layout separate \\
        --write-tables ${args} 2>&1 | tee ${prefix}/reconplot.log
    mv ${prefix}/per_chromosome/*.reconplot_{cn,sv}.tsv ${prefix}/

    Rscript ${reconplot_src}/run_reconplot.R ${source_args} --sample ${sample} --prefix ${sample} \\
        --genome ${genome} --outdir ${prefix}/genome_wide --regions all --layout together \\
        ${args} ${args2} 2>&1 | tee -a ${prefix}/reconplot.log
    ${focus_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reconplot: \$(Rscript -e 'cat(as.character(packageVersion("ReConPlot")))' 2>/dev/null)
        reconplot_wrapper: ${params.reconplot_dir ? 'local checkout' : params.reconplot_url}
        r-base: \$(Rscript -e 'cat(R.version\$major, R.version\$minor, sep=".")' 2>/dev/null)
        ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: (cn_source == sv_source ? "${cn_source}" : "${cn_source}_${sv_source}")
    """
    mkdir -p ${prefix}/per_chromosome ${prefix}/genome_wide
    touch ${prefix}/per_chromosome/${meta.id}_chr1.pdf ${prefix}/per_chromosome/${meta.id}_chr1.png
    touch ${prefix}/genome_wide/${meta.id}_genome_wide.pdf ${prefix}/genome_wide/${meta.id}_genome_wide.png
    touch ${prefix}/${meta.id}.reconplot_cn.tsv ${prefix}/${meta.id}.reconplot_sv.tsv ${prefix}/reconplot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reconplot: stub
        reconplot_wrapper: stub
        r-base: stub
        ggplot2: stub
    END_VERSIONS
    """
}
