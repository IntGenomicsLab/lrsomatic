process BCFTOOLS_NORM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'quay.io/biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf.gz"
    """
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --output ${prefix}.${extension} \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf}
    """

    stub:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf.gz"
    def index = ''
    if (extension in ['vcf.gz', 'bcf', 'bcf.gz']) {
        if (['--write-index=tbi', '-W=tbi'].any { arg -> args.contains(arg) } && extension == 'vcf.gz') {
            index = 'tbi'
        }
        else if (['--write-index=tbi', '-W=tbi', '--write-index=csi', '-W=csi', '--write-index', '-W'].any { arg -> args.contains(arg) }) {
            index = 'csi'
        }
    }
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = index ? "touch ${prefix}.${extension}.${index}" : ""
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
