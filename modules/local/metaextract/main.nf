process METAEXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), env(basecall_model),env(kinetics)  , emit: meta_ext
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ont = meta.platform == 'ont'
    """
    basecall_model=""
    kinetics=""
    if [ $ont = 'true' ]; then
        basecall_model=\$(samtools view -H ${bam} ${args} | awk -F'basecall_model=' '{print \$2}' | awk '{print \$1}'| tr -d '[:space:]')
    else
        kinetics=\$(samtools view -H ${bam} | awk '/--keep-kinetics/ {found=1} END {print (found ? "true" : "false")}')
        basecall_model="hifi_revio"
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
