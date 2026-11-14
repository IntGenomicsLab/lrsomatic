process DEEPSOMATIC_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: index
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_deepsomatic_filtered"
    """
    # DeepSomatic 1.10 declares FORMAT/NAF as Number=R but writes one value per ALT allele, so bcftools
    # cannot join two records at one position (norm -m +any). Declare it Number=A when the records say so.
    bcftools view -h ${vcf} > header.txt
    naf=\$( { bcftools query -i 'FMT/NAF!="."' -f '%ALT\\t[%NAF]\\n' ${vcf} || true; } | head -n 1 )
    if [ -n "\${naf}" ] && grep -q '^##FORMAT=<ID=NAF,Number=R,' header.txt; then
        n_alt=\$(echo "\${naf}" | cut -f 1 | awk -F, '{ print NF }')
        n_naf=\$(echo "\${naf}" | cut -f 2 | awk -F, '{ print NF }')
        if [ "\${n_naf}" -eq "\${n_alt}" ]; then
            sed -i 's/^##FORMAT=<ID=NAF,Number=R,/##FORMAT=<ID=NAF,Number=A,/' header.txt
        fi
    fi

    bcftools reheader -h header.txt ${vcf} \\
        | bcftools view \\
            --no-version \\
            ${args} \\
            -Oz \\
            --write-index=tbi \\
            -o ${prefix}.vcf.gz \\
            -
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_deepsomatic_filtered"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
