
process VCFSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(snv_vcf), path(indel_vcf)

    output:
    tuple val(meta), path("*somatic.vcf.gz"), emit: somatic_vcf
    tuple val(meta), path("*somatic.vcf.gz.tbi"), emit: somatic_tbi
    tuple val(meta), path("*germline.vcf.gz"), emit: germline_vcf
    tuple val(meta), path("*germline.vcf.gz.tbi"), emit: germline_tbi

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    bcftools view -i 'FILTER="PASS"' $indel_vcf | bgzip -c > indels_pass.vcf.gz
    bcftools view -i 'FILTER="PASS"' $snv_vcf | bgzip -c > snv_pass.vcf.gz
    tabix -p vcf indels_pass.vcf.gz
    tabix -p vcf snv_pass.vcf.gz
    bcftools concat -a -Oz -o somatic.vcf.gz indels_pass.vcf.gz snv_pass.vcf.gz
    tabix -p vcf somatic.vcf.gz

    bcftools view -i 'FILTER="NonSomatic"' $indel_vcf | bgzip -c > indels_filtered.vcf.gz
    bcftools view -i 'FILTER="NonSomatic"' $snv_vcf | bgzip -c > snv_filtered.vcf.gz
    tabix -p vcf indels_filtered.vcf.gz
    tabix -p vcf snv_filtered.vcf.gz
    bcftools concat -a -Oz -o germline_tmp.vcf.gz indels_filtered.vcf.gz snv_filtered.vcf.gz
    tabix -p vcf germline_tmp.vcf.gz

    bcftools view germline_tmp.vcf.gz | awk 'BEGIN{FS=OFS="\t"} /^#/ {print} !/^#/ { \$7="PASS"; print }' | \
        bgzip -c > germline.vcf.gz
    tabix -p vcf germline.vcf.gz

    # Cleanup intermediate files
    rm indels_pass.vcf.gz snv_pass.vcf.gz
    rm indels_pass.vcf.gz.tbi snv_pass.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfsplit: \$(bcftools --version |& sed '1!d ; s/bcftools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > somatic.vcf.gz
    echo "" | gzip > germline.vcf.gz
    echo "" | gzip > somatic.vcf.gz.tbi
    echo "" | gzip > germline.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfsplit: \$(bcftools --version |& sed '1!d ; s/bcftools //')
    END_VERSIONS
    """
}
