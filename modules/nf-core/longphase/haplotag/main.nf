process LONGPHASE_HAPLOTAG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83fce1d397cf71705cc096fc0e0e52f7013bdd471ef68ee53ae765688e5c439c/data':
        'community.wave.seqera.io/library/longphase_samtools:8c61296cae7a5fc0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snps), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.{bam,cram}"), emit: bam
    tuple val(meta), path("*.log")       , emit: log , optional: true
    tuple val("${task.process}"), val("longphase"), eval("longphase --version | head -n 1 | sed 's/Version: //'"), emit: versions_longphase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sv_file = svs ? "--sv-file ${svs}" : ""
    def mod_file = mods ? "--mod-file ${mods}" : ""
    // ext.cram: longphase writes its BAM into a FIFO that samtools encodes as CRAM (options from ext.args2,
    // e.g. embed_ref=1, which longphase's own --cram cannot set); no intermediate BAM touches the disk
    def cram_out = task.ext.cram ?: false
    def cram_start = cram_out ? """
    mkfifo ${prefix}.bam
    samtools view --threads ${task.cpus} -C -T ${fasta} ${args2} -o ${prefix}.cram ${prefix}.bam &
    conv=\$!
    trap 'kill \$conv 2>/dev/null || true' EXIT
    """ : ''
    def cram_end = cram_out ? """
    wait \$conv
    trap - EXIT
    rm ${prefix}.bam
    """ : ''

    """
    ${cram_start}
    longphase \\
        haplotag \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        --snp-file ${snps} \\
        --bam ${bam} \\
        ${sv_file} \\
        ${mod_file}
    ${cram_end}
    if [ -f "${prefix}.out" ]; then
        mv ${prefix}.out ${prefix}.log
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains('--cram') || task.ext.cram ? "cram" : "bam"
    def log = args.contains('--log') ? "touch ${prefix}.log" : ''
    """
    touch ${prefix}.${suffix}
    ${log}
    """
}
