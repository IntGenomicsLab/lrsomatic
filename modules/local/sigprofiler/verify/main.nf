process SIGPROFILER_VERIFY {
    tag "$genome"
    label 'process_single'

    // No conda: the image uses CHM13-T2T forks of SigProfilerMatrixGenerator (#250) and SigProfilerAssignment; see meta.yml
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-7894689'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-7894689'}"

    input:
    path(volume, stageAs: 'genome_volume')  // SigProfilerMatrixGenerator volume containing tsb/<genome>/
    val(genome)                             // SigProfilerMatrixGenerator genome name, e.g. GRCh38 or CHM13-T2T

    output:
    val(true)                                                                                                                , emit: verified
    tuple val("${task.process}"), val('sigprofilermatrixgenerator'), eval("python -c 'import importlib.metadata as m; print(m.version(\"SigProfilerMatrixGenerator\"))'"), topic: versions, emit: versions_sigprofilermatrixgenerator

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SIGPROFILER_VERIFY does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    """
    # SIGPROFILER_MATRIXGENERATOR re-checks the payload for every sample; checking once here fails a stale or damaged
    # volume before any sample work, with the reinstall instructions instead of a per-sample checksum error
    python - <<'PY'
    import sys
    from SigProfilerMatrixGenerator.scripts import reference_genome_manager as rgm

    manager = rgm.ReferenceGenomeManager("genome_volume")
    if not manager.is_genome_installed("${genome}"):
        manager.print_genome_checksum_verification_report("${genome}")
        sys.exit(
            "ERROR: the ${genome} payload in --sigprofiler_genome_dir does not match the checksums of this pipeline's "
            "SigProfilerMatrixGenerator. GRCh38 and CHM13-T2T payloads installed before lrsomatic PR #216 are a "
            "superseded revision. Reinstall it with --download_sigprofiler_genome (published to "
            "<outdir>/cache/sigprofiler/volume) and pass that directory on later runs."
        )
    PY
    """

    stub:
    """
    echo "stub: skipping checksum verification of genome_volume/tsb/${genome}"
    """
}
