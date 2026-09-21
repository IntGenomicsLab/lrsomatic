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
    tuple val(meta), path("*somatic.vcf.gz")        , emit: somatic_vcf
    tuple val(meta), path("*somatic.vcf.gz.tbi")    , emit: somatic_tbi
    tuple val(meta), path("*germline.vcf.gz")       , emit: germline_vcf
    tuple val(meta), path("*germline.vcf.gz.tbi")   , emit: germline_tbi

    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d ; s/bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    bcftools view -i 'FILTER="PASS"' $indel_vcf | bgzip -c > indels_pass.vcf.gz
    bcftools view -i 'FILTER="PASS"' $snv_vcf | bgzip -c > snv_pass.vcf.gz
    tabix -p vcf indels_pass.vcf.gz
    tabix -p vcf snv_pass.vcf.gz
    bcftools concat -a -Oz -o somatic_tmp.vcf.gz indels_pass.vcf.gz snv_pass.vcf.gz
    tabix -p vcf somatic_tmp.vcf.gz

    # Record the caller's original FILTER in INFO. These records are already PASS, but stamping
    # both splits keeps them symmetric and self-describing alongside the germline arm below.
    bcftools view somatic_tmp.vcf.gz | awk 'BEGIN{FS=OFS="\t"}
        /^##/ { print; next }
        /^#CHROM/ { print "##INFO=<ID=ORIG_FILTER,Number=1,Type=String,Description=\"FILTER value in the ClairS-TO output\">"; print; next }
        { \$8 = (\$8 == "." || \$8 == "") ? "ORIG_FILTER=" \$7 : \$8 ";ORIG_FILTER=" \$7; print }
    ' | bgzip -c > somatic.vcf.gz
    tabix -p vcf somatic.vcf.gz

    bcftools view -i 'FILTER~"NonSomatic" || INFO/Verdict_Germline=1' $indel_vcf | bgzip -c > indels_filtered.vcf.gz
    bcftools view -i 'FILTER~"NonSomatic" || INFO/Verdict_Germline=1' $snv_vcf | bgzip -c > snv_filtered.vcf.gz
    tabix -p vcf indels_filtered.vcf.gz
    tabix -p vcf snv_filtered.vcf.gz
    bcftools concat -a -Oz -o germline_tmp.vcf.gz indels_filtered.vcf.gz snv_filtered.vcf.gz
    tabix -p vcf germline_tmp.vcf.gz

    # FILTER is normalised to PASS so that downstream tools which filter on PASS -- implicitly or
    # otherwise -- see every germline record. The caller's original FILTER (typically NonSomatic)
    # would otherwise be destroyed here, which is what previously made germline records
    # indistinguishable from somatic ones once the two sets were merged for phasing; it is kept in
    # INFO/ORIG_FILTER instead. Germline/somatic provenance itself is stamped later, by
    # PHASING_HAPLOTYPING:TAG_GERMLINE / TAG_SOMATIC, which covers callers that bypass VCFSPLIT.
    bcftools view germline_tmp.vcf.gz | awk 'BEGIN{FS=OFS="\t"}
        /^##/ { print; next }
        /^#CHROM/ { print "##INFO=<ID=ORIG_FILTER,Number=1,Type=String,Description=\"FILTER value in the ClairS-TO output\">"; print; next }
        { \$8 = (\$8 == "." || \$8 == "") ? "ORIG_FILTER=" \$7 : \$8 ";ORIG_FILTER=" \$7; \$7 = "PASS"; print }
    ' | bgzip -c > germline.vcf.gz
    tabix -p vcf germline.vcf.gz

    # Cleanup intermediate files
    rm indels_pass.vcf.gz snv_pass.vcf.gz somatic_tmp.vcf.gz
    rm indels_pass.vcf.gz.tbi snv_pass.vcf.gz.tbi somatic_tmp.vcf.gz.tbi
    """

    stub:
    """
    echo "" | gzip > somatic.vcf.gz
    echo "" | gzip > germline.vcf.gz
    echo "" | gzip > somatic.vcf.gz.tbi
    echo "" | gzip > germline.vcf.gz.tbi
    """
}
