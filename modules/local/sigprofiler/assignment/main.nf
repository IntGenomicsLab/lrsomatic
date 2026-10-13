process SIGPROFILER_ASSIGNMENT {
    tag "$meta.id"
    label 'process_low'

    // No conda: the image uses CHM13-T2T forks of SigProfilerMatrixGenerator (#250) and SigProfilerAssignment; see meta.yml
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-7894689'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-7894689'}"

    input:
    tuple val(meta), path(vcf), path(sbs96), path(dbs78), path(id83)   // VCF to fit, plus its SigProfilerMatrixGenerator matrices (any may be []) deciding which contexts are fitted
    tuple val(meta2), path(volume, stageAs: 'genome_volume')           // SigProfilerMatrixGenerator volume containing tsb/<genome>/
    val(genome)                                                        // genome build, e.g. GRCh38 or CHM13-T2T
    val(cosmic_version)                                                // COSMIC version, e.g. 3.6

    output:
    tuple val(meta), path("*/SBS96")                                                                              , emit: sbs96                     , optional: true
    tuple val(meta), path("*/DBS78")                                                                              , emit: dbs78                     , optional: true
    tuple val(meta), path("*/ID83")                                                                               , emit: id83                      , optional: true
    tuple val(meta), path("*/*/Assignment_Solution/Activities/*_Activities.txt")                                  , emit: activities                , optional: true
    tuple val(meta), path("*/*/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/*.txt")           , emit: probabilities_per_mutation, optional: true
    tuple val("${task.process}"), val('sigprofilerassignment'), eval("python -c 'import importlib.metadata as m; print(m.version(\"SigProfilerAssignment\"))'"), topic: versions, emit: versions_sigprofilerassignment

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SIGPROFILER_ASSIGNMENT does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // SigProfilerMatrixGenerator names the sample from the file name up to the first '.'
    def sample = meta.id.toString().replaceAll(/[^A-Za-z0-9_-]/, '_')
    def reader = vcf.name.endsWith('.gz') ? 'gzip -cd' : 'cat'
    // VCF input is required for per-mutation probabilities; each context gets its own input dir for the internal matrix run
    // collapse_to_SBS96 defaults to True, which SigProfilerAssignment refuses for DBS78 and ID83
    def fits = [[sbs96, 'SBS96', '96', 'True'], [dbs78, 'DBS78', 'DINUC', 'False'], [id83, 'ID83', 'ID', 'False']]
        .findAll { matrix, _ctx, _context_type, _collapse -> matrix }
        .collect { matrix, ctx, context_type, collapse ->
            """
    if [ "\$(matrix_total ${matrix})" -gt 0 ]; then
        mkdir -p input_${ctx}
        ${reader} ${vcf} > input_${ctx}/${sample}.vcf
        SigProfilerAssignment cosmic_fit \\
            input_${ctx}/ \\
            ${prefix}/${ctx} \\
            --input_type vcf \\
            --context_type ${context_type} \\
            --collapse_to_SBS96 ${collapse} \\
            --genome_build ${genome} \\
            --cosmic_version ${cosmic_version} \\
            --export_probabilities_per_mutation True \\
            --cpu ${task.cpus} \\
            ${args}
    else
        echo "${ctx}: no mutations in ${matrix}, skipping COSMIC fit" >&2
    fi
    """
        }
        .join('\n')
    """
    # total mutation count of a SigProfilerMatrixGenerator matrix (all samples, all classes)
    matrix_total() {
        awk -F'\\t' 'NR > 1 { for (i = 2; i <= NF; i++) s += \$i } END { printf "%d\\n", s }' "\$1"
    }

    # Symlink the genome payload into a local volume so the tool never writes into the shared reference dir
    mkdir -p ${prefix} volume/tsb plot_templates
    ln -s "\$(readlink -f ${volume})/tsb/${genome}" volume/tsb/${genome}
    # cosmic_fit ignores --volume, so the genome and template dirs are set through the environment
    export SIGPROFILERMATRIXGENERATOR_VOLUME="\$PWD/volume"
    export SIGPROFILERASSIGNMENT_VOLUME="\$PWD/volume"
    # sigProfilerPlotting caches templates in its read-only package dir unless this points elsewhere
    export SIGPROFILERPLOTTING_VOLUME="\$PWD/plot_templates"
    ${fits}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample = meta.id.toString().replaceAll(/[^A-Za-z0-9_-]/, '_')
    """
    mkdir -p ${prefix}/SBS96/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities
    touch ${prefix}/SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt
    touch ${prefix}/SBS96/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_${sample}.txt
    """
}
