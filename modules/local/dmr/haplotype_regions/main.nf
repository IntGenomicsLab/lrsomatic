process DMR_HAPLOTYPE_REGIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    // hp1/hp2 bedMethyl are not read for methylation values here, only for which CpG positions
    // they cover -- restricting cpg_islands_bed to islands covered in BOTH haplotypes is what
    // keeps modkit dmr pair from comparing real data on one haplotype against the other
    // haplotype's absence of data at the same island.
    tuple val(meta), path(hp1_bedmethyl), path(hp1_tbi), path(hp2_bedmethyl), path(hp2_tbi)
    path cpg_islands_bed
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.regions.bed"), emit: regions_bed
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f1,2 ${fai} > genome.txt

    zcat ${hp1_bedmethyl} | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}' | sort -k1,1 -k2,2n -u > hp1.cpg.bed
    zcat ${hp2_bedmethyl} | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}' | sort -k1,1 -k2,2n -u > hp2.cpg.bed

    bedtools intersect -u -a ${cpg_islands_bed} -b hp1.cpg.bed | sort -k1,1 -k2,2n > islands.hp1.bed
    bedtools intersect -u -a ${cpg_islands_bed} -b hp2.cpg.bed | sort -k1,1 -k2,2n > islands.hp2.bed
    bedtools intersect -u -a islands.hp1.bed -b islands.hp2.bed | bedtools sort -g genome.txt -i - > ${prefix}.regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
