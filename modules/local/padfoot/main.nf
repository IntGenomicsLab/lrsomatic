process PADFOOT {
    tag "${meta.id}:${sv_caller}+${cna_caller}"
    label 'process_medium'

    // No conda: the image ships Padfoot itself (a pinned Tim-Yu/Padfoot commit at /opt/padfoot), not just its
    // dependencies (guard in `script:`). Built from containers/padfoot/Dockerfile: Padfoot + RepeatMasker 4.2.4 with the
    // Dfam 4.0 root and curated-consensus partitions. Padfoot update = new commit in the Dockerfile, rebuild, bump these two tags.
    // Override per site with `process { withName: '.*:PADFOOT_(SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/timmy9527/padfoot-repeatmasker-sif:4.2.4-dfam4-padfoot-1748e84-r2'
        : 'docker.io/timmy9527/padfoot-repeatmasker:4.2.4-dfam4-padfoot-1748e84-r2'}"

    input:
    // ploidy_file: the CN caller's fitted purity/ploidy table (SAVANA *_fitted_purity_ploidy.tsv, Wakhan
    // solutions_ranks.tsv) or []; Padfoot labels gene copy number against ploidy/2 and otherwise estimates
    // the ploidy from the profile with a warning.
    tuple val(meta), path(sv_vcf), val(sv_caller), path(cna_file), val(cna_caller), path(ploidy_file)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), val(genome), path(gff), path(rm) // gff/rm may be [] -> Padfoot bundled annotations for `genome`

    output:
    tuple val(meta), path("${prefix}/annotated_svs.tsv"), emit: annotated_svs
    tuple val(meta), path("${prefix}/by_gene.tsv")      , emit: by_gene
    tuple val(meta), path("${prefix}/padfoot.log")      , emit: log
    tuple val(meta), path("${prefix}/staged.txt")       , emit: staged, optional: true // stub only: the staged inputs, for the tests
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PADFOOT does not support Conda: Padfoot ships only inside its container. Use Docker / Singularity / Apptainer, or --skip_padfoot."
    }
    def padfoot = '/opt/padfoot'   // Padfoot source tree inside the image (padfoot.py + beds/)
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${sv_caller}_${cna_caller}"
    def gff_arg = gff ? "--gff ${gff}" : ''
    def rm_arg  = rm  ? "--rm ${rm}"   : ''
    def ploidy_arg = ploidy_file ? "--ploidy-file ${ploidy_file}" : ''

    """
    python3 ${padfoot}/padfoot.py \\
        --sv-vcf ${sv_vcf} \\
        --sv-caller ${sv_caller} \\
        --cna-file ${cna_file} \\
        --cna-caller ${cna_caller} \\
        --ref ${fasta} \\
        --genome ${genome} \\
        ${gff_arg} \\
        ${rm_arg} \\
        ${ploidy_arg} \\
        --threads ${task.cpus} \\
        --out-dir ${prefix} \\
        ${args}

    rm -rf ${prefix}/temp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        padfoot: \$(python3 ${padfoot}/padfoot.py --version 2>&1 | tail -1)
        padfoot_commit: \${PADFOOT_COMMIT:-unknown}
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        repeatmasker: \$(command -v RepeatMasker >/dev/null && RepeatMasker -v 2>&1 | sed -n 's/^RepeatMasker version //p' || echo 'not available')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${sv_caller}_${cna_caller}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/annotated_svs.tsv ${prefix}/by_gene.tsv ${prefix}/padfoot.log
    printf '%s\\n' ${sv_vcf} ${cna_file} ${ploidy_file} > ${prefix}/staged.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        padfoot: stub
        padfoot_commit: stub
        minimap2: stub
        samtools: stub
        repeatmasker: stub
    END_VERSIONS
    """
}
