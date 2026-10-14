process SHATTERSEEK {
    tag "${meta.id}:${cn_source}+${sv_source}"
    label 'process_low'

    // No conda: the ShatterSeek fork (chm13, corrected statistics, PCAWG tiers) ships only in this image, built from
    // containers/shatterseek/Dockerfile; the wrapper is pipeline glue in assets/shatterseek/ and reads the callers with
    // the ReConPlot parsers in assets/reconplot/. Fork update = new commit there, rebuild, bump this tag.
    container "docker.io/robertaforsyth/shatterseek:1.1.1-chm13-e1cee6c"

    input:
    // cn_files / sv_files are the caller output files the parsers discover by name.
    // When cn_source == sv_source (savana) put everything in cn_files and leave sv_files empty.
    tuple val(meta), val(cn_source), path(cn_files, stageAs: 'cn_input/*'), val(sv_source), path(sv_files, stageAs: 'sv_input/*')
    tuple val(meta2), path(shatterseek_src)   // ShatterSeek wrapper (run_shatterseek.R + VERSION)
    tuple val(meta3), path(reconplot_src)     // ReConPlot wrapper, for its caller parsers (R/)
    val(genome)                               // hg38 | chm13 | hg19

    output:
    tuple val(meta), path("${prefix}/${meta.id}.shatterseek.tsv")               , emit: summary
    tuple val(meta), path("${prefix}/${meta.id}.shatterseek_{sv,cn}_input.tsv") , emit: inputs
    tuple val(meta), path("${prefix}/plots/*.pdf")                              , emit: plots, optional: true
    tuple val(meta), path("${prefix}/shatterseek.log")                          , emit: log
    tuple val(meta), path("${prefix}/staged.txt")                               , emit: staged, optional: true // stub only: the staged inputs, for the tests
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SHATTERSEEK does not support Conda: the ShatterSeek fork ships only inside its container. Use Docker / Singularity / Apptainer, or --skip_shatterseek."
    }
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: (cn_source == sv_source ? "${cn_source}" : "${sv_source}_${cn_source}")
    def sv_in = cn_source == sv_source ? '' : '--sv-input sv_input'
    """
    mkdir -p ${prefix}

    Rscript ${shatterseek_src}/run_shatterseek.R \\
        --reconplot-dir ${reconplot_src} \\
        --cn-source ${cn_source} --sv-source ${sv_source} --cn-input cn_input ${sv_in} \\
        --sample ${meta.id} --prefix ${meta.id} --genome ${genome} --outdir ${prefix} \\
        ${args} 2>&1 | tee ${prefix}/shatterseek.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shatterseek: \$(Rscript -e 'cat(as.character(packageVersion("ShatterSeek")))' 2>/dev/null)
        shatterseek_fork_commit: \${SHATTERSEEK_COMMIT:-unknown}
        shatterseek_wrapper: \$(cat ${shatterseek_src}/VERSION 2>/dev/null || echo unknown)
        r-base: \$(Rscript -e 'cat(R.version\$major, R.version\$minor, sep=".")' 2>/dev/null)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: (cn_source == sv_source ? "${cn_source}" : "${sv_source}_${cn_source}")
    """
    mkdir -p ${prefix}/plots
    touch ${prefix}/${meta.id}.shatterseek.tsv ${prefix}/${meta.id}.shatterseek_sv_input.tsv ${prefix}/${meta.id}.shatterseek_cn_input.tsv
    touch ${prefix}/plots/${meta.id}_chr1.pdf ${prefix}/shatterseek.log
    { ls -1 cn_input/* sv_input/* 2>/dev/null || true; } | sort > ${prefix}/staged.txt   # sv_input is absent in single-source mode

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shatterseek: stub
        shatterseek_fork_commit: stub
        shatterseek_wrapper: stub
        r-base: stub
    END_VERSIONS
    """
}
