process VCFSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'biocontainers/bcftools:1.24--h118bc1c_2' }"

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
    // Extra bcftools view filter for the germline arm (--smallvar_filter_pass sets it in conf/modules.config)
    def germline_args = task.ext.args ?: ''
    """

    bcftools view -i 'FILTER="PASS"' $indel_vcf | bgzip -c > indels_pass.vcf.gz
    bcftools view -i 'FILTER="PASS"' $snv_vcf | bgzip -c > snv_pass.vcf.gz
    tabix -p vcf indels_pass.vcf.gz
    tabix -p vcf snv_pass.vcf.gz
    bcftools concat -a -Oz -o somatic.vcf.gz indels_pass.vcf.gz snv_pass.vcf.gz
    tabix -p vcf somatic.vcf.gz

    bcftools view -i 'FILTER~"NonSomatic" || INFO/Verdict_Germline=1' $indel_vcf | bcftools view ${germline_args} - | bgzip -c > indels_filtered.vcf.gz
    bcftools view -i 'FILTER~"NonSomatic" || INFO/Verdict_Germline=1' $snv_vcf | bcftools view ${germline_args} - | bgzip -c > snv_filtered.vcf.gz
    tabix -p vcf indels_filtered.vcf.gz
    tabix -p vcf snv_filtered.vcf.gz
    bcftools concat -a -Oz -o germline_tmp.vcf.gz indels_filtered.vcf.gz snv_filtered.vcf.gz
    tabix -p vcf germline_tmp.vcf.gz

    # Normalise FILTER to PASS, keeping the original in INFO/ORIG_FILTER (";" stored as ",").
    # An ORIG_FILTER already present (header or record) is kept rather than duplicated.
    bcftools view germline_tmp.vcf.gz | awk -v q='"' 'BEGIN{FS=OFS="\t"}
        /^##INFO=<ID=ORIG_FILTER,/ { has_hdr = 1 }
        /^##/ { print; next }
        /^#CHROM/ { if (!has_hdr) print "##INFO=<ID=ORIG_FILTER,Number=.,Type=String,Description=" q "Original FILTER value before normalisation to PASS" q ">"; print; next }
        { of = \$7; gsub(/;/, ",", of)
          if (\$8 == "." || \$8 == "") \$8 = "ORIG_FILTER=" of
          else if (\$8 !~ /(^|;)ORIG_FILTER=/) \$8 = \$8 ";ORIG_FILTER=" of
          \$7 = "PASS"
          print }
    ' | bgzip -c > germline.vcf.gz
    tabix -p vcf germline.vcf.gz

    # Fail here, not downstream, if either header does not parse.
    bcftools view -h somatic.vcf.gz > /dev/null
    bcftools view -h germline.vcf.gz > /dev/null

    # Cleanup intermediate files
    rm indels_pass.vcf.gz snv_pass.vcf.gz
    rm indels_pass.vcf.gz.tbi snv_pass.vcf.gz.tbi
    """

    stub:
    """
    echo "" | gzip > somatic.vcf.gz
    echo "" | gzip > germline.vcf.gz
    echo "" | gzip > somatic.vcf.gz.tbi
    echo "" | gzip > germline.vcf.gz.tbi
    """
}
