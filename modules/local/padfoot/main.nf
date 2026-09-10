process PADFOOT {
    tag "${meta.id}:${sv_caller}+${cna_caller}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/padfoot/Dockerfile: Padfoot deps + RepeatMasker 4.2.4 + Dfam 4.0 (Padfoot itself is staged as source).
    // Override per site with `process { withName: '.*:PADFOOT_(SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.
    container "ghcr.io/tim-yu/padfoot-repeatmasker@sha256:f98b0d352ec47cd9fa015f321959dcba7f6c6cd4da05beaea8d811b97dff70f5"

    input:
    tuple val(meta), path(sv_vcf), val(sv_caller), path(cna_file), val(cna_caller)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(padfoot_src)              // Padfoot source tree (contains padfoot.py + beds/)
    tuple val(meta5), val(genome), path(gff), path(rm) // gff/rm may be [] -> Padfoot bundled annotations for `genome`

    output:
    tuple val(meta), path("${prefix}/annotated_svs.tsv"), emit: annotated_svs
    tuple val(meta), path("${prefix}/by_gene.tsv")      , emit: by_gene
    tuple val(meta), path("${prefix}/padfoot.log")      , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${sv_caller}_${cna_caller}"
    def gff_arg = gff ? "--gff ${gff}" : ''
    def rm_arg  = rm  ? "--rm ${rm}"   : ''

    """
    python3 ${padfoot_src}/padfoot.py \\
        --sv-vcf ${sv_vcf} \\
        --sv-caller ${sv_caller} \\
        --cna-file ${cna_file} \\
        --cna-caller ${cna_caller} \\
        --ref ${fasta} \\
        --genome ${genome} \\
        ${gff_arg} \\
        ${rm_arg} \\
        --threads ${task.cpus} \\
        --out-dir ${prefix} \\
        ${args}

    rm -rf ${prefix}/temp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        padfoot: \$(python3 ${padfoot_src}/padfoot.py --version 2>&1 | tail -1)
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        padfoot: stub
        minimap2: stub
        samtools: stub
        repeatmasker: stub
    END_VERSIONS
    """
}
