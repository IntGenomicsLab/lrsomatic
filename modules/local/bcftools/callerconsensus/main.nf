process BCFTOOLS_CALLER_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(prio_vcf), path(prio_tbi), path(other_vcf), path(other_tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_consensus"
    // Which of the other VCF's records (ext.args2, before splitting) and split alleles (ext.args, after) count as
    // matches; by default the alleles its genotypes carry. ext.args = '' matches every allele of the selected records.
    def args  = task.ext.args != null ? task.ext.args : "-i 'GT=\"alt\"'"
    def args2 = task.ext.args2 ?: ''
    """
    # The called alleles of the priority caller, one allele per record: split, left-align, keep alleles the genotype carries.
    # Its records are numbered first (in ID), so each split allele points back to its record.
    bcftools view --no-version ${prio_vcf} \\
        | awk 'BEGIN { FS = OFS = "\\t" } /^#/ { print; next } { \$3 = ++n; print }' \\
        | bcftools norm --no-version -f ${fasta} -m -any --check-ref w -Ou - \\
        | bcftools sort -T ./ -Ou \\
        | bcftools view --no-version -i 'GT="alt"' -Oz -W=tbi -o prio_called.vcf.gz
    # The other VCF's matching alleles, likewise split and left-aligned
    bcftools view --no-version ${args2} -Ou ${other_vcf} \\
        | bcftools norm --no-version -f ${fasta} -m -any --check-ref w -Ou - \\
        | bcftools sort -T ./ -Ou \\
        | bcftools view --no-version ${args} -Oz -W=tbi -o other_called.vcf.gz

    # Records with a called allele the other VCF matches (same CHROM, POS, REF, ALT), and those with one it does not
    bcftools isec --no-version -c none -n =2 -w 1 prio_called.vcf.gz other_called.vcf.gz -Ou \\
        | bcftools query -f '%ID\\n' > shared.ids
    bcftools isec --no-version -c none -n =1 -w 1 prio_called.vcf.gz other_called.vcf.gz -Ou \\
        | bcftools query -f '%ID\\n' > prio_only.ids

    # Output the priority caller's own records, unchanged, where every called allele is matched
    bcftools view --no-version ${prio_vcf} \\
        | awk 'BEGIN { FS = OFS = "\\t" }
            FILENAME == "shared.ids"    { shared[\$1] = 1; next }
            FILENAME == "prio_only.ids" { prio_only[\$1] = 1; next }
            /^#/ { print; next }
            { n++; if ((n in shared) && !(n in prio_only)) print }' shared.ids prio_only.ids - \\
        | bcftools view --no-version -Oz -W=tbi -o ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_consensus"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
