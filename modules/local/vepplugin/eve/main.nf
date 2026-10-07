process VEPPLUGIN_EVE {
    tag "${eve_dir}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    path eve_dir

    output:
    path "eve_merged.vcf.gz{,.tbi}", emit: files
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # TMPDIR into the task dir: it is on scratch, and an inherited host TMPDIR need not be bound in the container
    export TMPDIR=\$PWD

    prepare_vep_plugin_data.sh eve ${eve_dir} . ${args}
    """

    stub:
    """
    touch eve_merged.vcf.gz
    touch eve_merged.vcf.gz.tbi
    """
}
