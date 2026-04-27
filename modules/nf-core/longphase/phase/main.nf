process LONGPHASE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83fce1d397cf71705cc096fc0e0e52f7013bdd471ef68ee53ae765688e5c439c/data':
        'community.wave.seqera.io/library/longphase_samtools:8c61296cae7a5fc0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snvs), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("${prefix}.vcf.gz")        , emit: snv_vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")    , emit: snv_vcf_index
    tuple val(meta), path("${prefix}_SV.vcf.gz")     , emit: sv_vcf , optional: true
    tuple val(meta), path("${prefix}_SV.vcf.gz.tbi") , emit: sv_vcf_index , optional: true
    tuple val(meta), path("${prefix}_mod.vcf.gz")    , emit: mod_vcf, optional: true
    tuple val(meta), path("${prefix}_mod.vcf.gz.tbi"), emit: mod_vcf_index, optional: true
    tuple val("${task.process}"), val("longphase"), eval("longphase --version | head -n 1 | sed 's/Version: //'"), emit: versions_longphase, topic: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sv_file = svs ? "--sv-file ${svs}" : ""
    def mod_file = mods ? "--mod-file ${mods}" : ""
    def bams = bam.collectMany { file -> ["-b", file] }.join(" ")
    """
    longphase \\
        phase \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        --snp-file ${snvs} \\
        ${bams} \\
        ${sv_file} \\
        ${mod_file} \\

    bgzip \\
        --threads $task.cpus \\
        $args2 \\
        ${prefix}*.vcf

    tabix -p vcf ${prefix}.vcf.gz

    if [ -f ${prefix}_SV.vcf.gz ]; then
        tabix -p vcf ${prefix}_SV.vcf.gz
    fi

    if [ -f ${prefix}_mod.vcf.gz ]; then
        tabix -p vcf ${prefix}_mod.vcf.gz
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
      tabix -p vcf ${prefix}.vcf.gz

    if [ -f ${prefix}_SV.vcf.gz ]; then
        tabix -p vcf ${prefix}_SV.vcf.gz
    fi

    if [ -f ${prefix}_mod.vcf.gz ]; then
        tabix -p vcf ${prefix}_mod.vcf.gz
    fi
    """
}
