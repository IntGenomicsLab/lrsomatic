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

    genome = "${genome}"
    if genome not in rgm.CHECKSUMS:
        sys.exit(
            f"ERROR: this pipeline's SigProfilerMatrixGenerator has no checksums for {genome} (registered: "
            f"{', '.join(sorted(rgm.CHECKSUMS))}). Set --sigprofiler_genome to one of them or use --skip_signatures."
        )

    manager = rgm.ReferenceGenomeManager("genome_volume")
    if not manager.is_genome_installed(genome):
        manager.print_genome_checksum_verification_report(genome)
        tsb = manager.reference_dir.get_tsb_dir() / genome
        expected = rgm.CHECKSUMS[genome]
        missing = [chrom for chrom in expected if not (tsb / f"{chrom}.txt").is_file()]
        if missing:
            cause = (
                f"is an incomplete install: {len(missing)} of {len(expected)} chromosome files are missing "
                f"({', '.join(missing)})."
            )
        else:
            cause = (
                "does not match the checksums of this pipeline's SigProfilerMatrixGenerator. Either it is a stale "
                "payload (GRCh38 and CHM13-T2T installed for lrsomatic < 1.2.0 are a superseded revision) or the "
                "copy is corrupted."
            )
        sys.exit(
            f"ERROR: the {genome} payload in --sigprofiler_genome_dir {cause} Reinstall it with "
            "--download_sigprofiler_genome (published to <outdir>/cache/sigprofiler/volume) and pass that directory "
            "on later runs."
        )
    PY
    """

    stub:
    """
    echo "stub: skipping checksum verification of genome_volume/tsb/${genome}"
    """
}
