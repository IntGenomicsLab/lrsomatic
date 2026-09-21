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
    # The header line's double quotes arrive via -v q. Escaped quotes inside a Nextflow script
    # block are fragile, and losing the escape silently produces an unparseable VCF header.
    # FILTER is ";"-delimited but ";" separates INFO fields, so it is stored as ",".
    bcftools view somatic_tmp.vcf.gz | awk -v q='"' 'BEGIN{FS=OFS="\t"}
        /^##/ { print; next }
        /^#CHROM/ { print "##INFO=<ID=ORIG_FILTER,Number=.,Type=String,Description=" q "FILTER value in the ClairS-TO output" q ">"; print; next }
        { of = \$7; gsub(/;/, ",", of)
          \$8 = (\$8 == "." || \$8 == "") ? "ORIG_FILTER=" of : \$8 ";ORIG_FILTER=" of
          print }
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
    bcftools view germline_tmp.vcf.gz | awk -v q='"' 'BEGIN{FS=OFS="\t"}
        /^##/ { print; next }
        /^#CHROM/ { print "##INFO=<ID=ORIG_FILTER,Number=.,Type=String,Description=" q "FILTER value in the ClairS-TO output" q ">"; print; next }
        { of = \$7; gsub(/;/, ",", of)
          \$8 = (\$8 == "." || \$8 == "") ? "ORIG_FILTER=" of : \$8 ";ORIG_FILTER=" of
          \$7 = "PASS"
          print }
    ' | bgzip -c > germline.vcf.gz
    tabix -p vcf germline.vcf.gz

    # Read both headers back. tabix will happily index a VCF whose header htslib cannot parse, so
    # without this a malformed header surfaces as a confusing failure in a later process instead
    # of here. set -e is in effect, so a bad header fails this task.
    bcftools view -h somatic.vcf.gz > /dev/null
    bcftools view -h germline.vcf.gz > /dev/null

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
