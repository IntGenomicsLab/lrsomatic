process VCFTAG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val  flag

    output:
    tuple val(meta), path("${prefix}.vcf.gz")     , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi") , emit: tbi

    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d ; s/bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_${flag.toLowerCase()}"
    """
    # Stamp a constant INFO flag recording which arm this record came from, and normalise FILTER
    # to PASS. Provenance has to live in INFO rather than FILTER: downstream steps rewrite and
    # filter on FILTER, so a FILTER-based label does not survive to where it is needed.
    # bcftools annotate cannot set a constant INFO field without an annotation file, hence awk.
    bcftools view ${vcf} | awk -v flag="${flag}" -v q='"' 'BEGIN{FS=OFS="\t"}
        /^##/ { print; next }
        /^#CHROM/ {
            print "##INFO=<ID=" flag ",Number=0,Type=Flag,Description=" q "Record originates from the " tolower(flag) " call set" q ">"
            print
            next
        }
        {
            \$7 = "PASS"
            \$8 = (\$8 == "." || \$8 == "") ? flag : \$8 ";" flag
            print
        }
    ' | bgzip -c > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${flag.toLowerCase()}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
