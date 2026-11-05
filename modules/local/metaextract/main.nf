process METAEXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), env(basecall_model), env(kinetics)  , emit: meta_ext
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ont = meta.platform == 'ont'
    basecall_model = ''
    kinetics = ''
    """
    export basecall_model="${basecall_model}"
    export kinetics="${kinetics}"
    if [ $ont = 'true' ]; then
        basecall_model=\$(samtools view -H "${bam}" ${args} | awk -F'basecall_model=' '/basecall_model=/ {print \$2; exit}' | awk '{print \$1}' | tr -d '[:space:]')
    else
        kinetics=\$(samtools view -H ${bam} | awk '/--keep-kinetics/ {found=1} END {print (found ? "true" : "false")}')
        basecall_model="hifi_revio"
    fi
    """

    stub:
    // A real model name, so the Clair3/ClairS model lookups resolve and -stub runs reach the callers
    def stub_model = meta.platform == 'ont' ? 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0' : 'hifi_revio'
    def stub_kinetics = meta.platform == 'ont' ? '' : 'false'
    """
    export basecall_model="${stub_model}"
    export kinetics="${stub_kinetics}"
    """
}
