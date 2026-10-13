process SIGPROFILER_ASSIGNMENT {
    tag "$meta.id"
    label 'process_low'

    // No conda: the image uses CHM13-T2T forks of SigProfilerMatrixGenerator (#250) and SigProfilerAssignment; see meta.yml
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-7894689'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-7894689'}"

    input:
    tuple val(meta), path(spmg_output, stageAs: 'spmg/output'), path(sbs96), path(dbs78), path(id83)   // MatrixGenerator output dir and matrices (any may be [])
    val(genome)                                                                                         // genome build for the COSMIC reference signatures, e.g. GRCh38 or CHM13-T2T
    val(cosmic_version)                                                                                 // COSMIC version, e.g. 3.6

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
    // SBS matrices are collapsed to SBS96 for the fit; DBS78 and ID83 must not be
    // (SigProfilerAssignment refuses collapse_to_SBS96 for non-SBS contexts).
    def fits = [[sbs96, 'SBS96', 'True', '96'], [dbs78, 'DBS78', 'False', '78'], [id83, 'ID83', 'False', '83']]
        .findAll { matrix, _ctx, _collapse, _m -> matrix }
        .collect { matrix, ctx, collapse, m ->
            """
    if [ "\$(matrix_total ${matrix})" -gt 0 ]; then
        SigProfilerAssignment cosmic_fit \\
            ${matrix} \\
            ${prefix}/${ctx} \\
            --input_type matrix \\
            --genome_build ${genome} \\
            --cosmic_version ${cosmic_version} \\
            --collapse_to_SBS96 ${collapse} \\
            --volume spa_volume \\
            --cpu ${task.cpus} \\
            ${args}
        python per_mutation.py ${prefix}/${ctx} ${m}
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

    # cosmic_fit exports per-mutation probabilities only from VCF input, so join the fit with MatrixGenerator's seqinfo
    cat > per_mutation.py <<'PY'
    import os, sys
    import pandas as pd
    from SigProfilerAssignment.decompose_subroutines import probabilities_per_mutation
    fit, m = sys.argv[1:]
    activities = os.path.join(fit, "Assignment_Solution", "Activities")
    probs = pd.read_csv(os.path.join(activities, "Decomposed_MutationType_Probabilities.txt"), sep="\\t", index_col=0)
    out = os.path.join(activities, "Decomposed_Mutation_Probabilities")
    os.makedirs(out, exist_ok=True)
    # same index and sort as SigProfilerAssignment's own VCF-mode export
    for df, sample in zip(*probabilities_per_mutation(probs, "spmg", m)):
        df.set_index("Sample Names").sort_values(by=["Chr", "Pos"]).to_csv(os.path.join(out, f"Decomposed_Mutation_Probabilities_{sample}.txt"), sep="\\t")
    PY

    # sigProfilerPlotting caches templates in its read-only package dir unless this points elsewhere
    mkdir -p ${prefix} spa_volume
    export SIGPROFILERPLOTTING_VOLUME="\$PWD/spa_volume"
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
