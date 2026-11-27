process BCFTOOLS_CONCAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'quay.io/biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(vcfs), path(tbi)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: vcf
    tuple val(meta), path("${prefix}.${extension}.{tbi,csi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tbi_names = tbi.findAll { file -> !(file instanceof List) }.collect { file -> file.name }
    def create_input_index = vcfs.collect { vcf -> tbi_names.contains(vcf.name + ".tbi") || tbi_names.contains(vcf.name + ".csi") ? "" : "tabix ${vcf}" }.join("\n    ")
    extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def input = vcfs.sort { vcf -> vcf.toString() }.join(" ")
    """
    ${create_input_index}

    bcftools concat \\
        --output ${prefix}.${extension} \\
        ${args} \\
        --threads ${task.cpus} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index_extension = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index_extension.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index_extension}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
