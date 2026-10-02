process DMR_HAPLOTYPE_REGIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // The plain bedtools:2.31.1 biocontainer has no tabix binary at all -- needed since
    // haplotype_regions/main.nf queries the bgzip+tabix-indexed bedMethyl inputs directly.
    // Reusing nf-core/modules' pints/caller container here (real, already in production use,
    // confirmed via its own environment.yml to bundle both bedtools and htslib) rather than
    // hand-constructing an unverifiable Wave/mulled tag; it carries unrelated pybedtools/pypints
    // baggage this module doesn't use, but that's preferable to guessing a container hash.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a9e30012e1b41baf9acd1ff94e01161138d8aa17f4e97aa32f2dc4effafcd1/data'
        : 'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6'}"

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
    // Minimum valid coverage (bedMethyl column 10, Nvalid_cov) for a position to count as
    // "covered" -- 1 preserves the previous any-depth behaviour; override via ext.args to
    // require more.
    def min_valid_coverage = task.ext.args ?: 1
    """
    cut -f1,2 ${fai} > genome.txt

    # Query directly via the bgzip+tabix index for positions inside the CpG islands, instead
    # of decompressing/sorting every bedMethyl row genome-wide -- hp*_bedmethyl are already
    # coordinate-sorted and tabix-indexed upstream, so this only reads the (small) subset of
    # rows that can possibly matter here.
    tabix -R ${cpg_islands_bed} ${hp1_bedmethyl} \\
        | awk -v min_cov=${min_valid_coverage} 'BEGIN{OFS="\\t"} \$10>=min_cov {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n -u > hp1.cpg.bed
    tabix -R ${cpg_islands_bed} ${hp2_bedmethyl} \\
        | awk -v min_cov=${min_valid_coverage} 'BEGIN{OFS="\\t"} \$10>=min_cov {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n -u > hp2.cpg.bed

    bedtools intersect -u -a ${cpg_islands_bed} -b hp1.cpg.bed | sort -k1,1 -k2,2n > islands.hp1.bed
    bedtools intersect -u -a ${cpg_islands_bed} -b hp2.cpg.bed | sort -k1,1 -k2,2n > islands.hp2.bed
    bedtools intersect -u -a islands.hp1.bed -b islands.hp2.bed | bedtools sort -g genome.txt -i - > ${prefix}.regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//g')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //')
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
