process BCFTOOLS_CALLER_UNION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(prio_vcf), path(prio_tbi), path(other_vcf), path(other_tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_union"
    """
    # One record per CHROM:POS, per caller and across callers: PASS first, then the rest; the priority caller wins within each tier.
    bcftools view --no-version -f PASS ${prio_vcf} -Ou \\
        | bcftools norm --no-version -d all -Oz -W=tbi -o prio_pass.vcf.gz
    bcftools view --no-version -f PASS -T ^prio_pass.vcf.gz ${other_vcf} -Ou \\
        | bcftools norm --no-version -d all -Oz -W=tbi -o other_pass.vcf.gz
    bcftools concat --no-version -a prio_pass.vcf.gz other_pass.vcf.gz -Oz -W=tbi -o pass.vcf.gz

    bcftools view --no-version -e 'FILTER="PASS"' -T ^pass.vcf.gz ${prio_vcf} -Ou \\
        | bcftools norm --no-version -d all -Oz -W=tbi -o prio_rest.vcf.gz
    bcftools view --no-version -e 'FILTER="PASS"' -T ^pass.vcf.gz ${other_vcf} -Ou \\
        | bcftools view --no-version -T ^prio_rest.vcf.gz -Ou \\
        | bcftools norm --no-version -d all -Oz -W=tbi -o other_rest.vcf.gz

    bcftools concat --no-version -a pass.vcf.gz prio_rest.vcf.gz other_rest.vcf.gz -Ou \\
        | bcftools sort -T ./ -Oz -W=tbi -o ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_union"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
