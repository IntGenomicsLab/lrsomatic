process CLAIRSTO_CNA_RESOURCES {
    tag "$genome_name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04':
        'nf-core/ubuntu:22.04' }"

    input:
    path(loci_files)      // per-chromosome 1000G loci files for this reference
    path(allele_files)    // per-chromosome 1000G allele files for this reference
    path(gc_file)         // GC content per locus
    path(rt_file)         // replication timing per locus
    val(genome_name)      // reference build name, for the tag only

    output:
    path("clairsto_cna_resources"), emit: cna_resources
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ClairS-TO's Verdict module hard-codes the GRCh38 resource file names in
    # src/cna_germline_tagging.py (G1000_loci_hg38_<ctg>.txt, G1000_alleles_hg38_, GC_G1000_hg38.txt,
    # RT_G1000_hg38.txt) and has no reference-build switch. Present this reference's files under
    # those names so --cna_resource_dir picks them up. Verdict covers chr1-22 and chrX only.
    mkdir -p clairsto_cna_resources/loci_files clairsto_cna_resources/allele_files

    for ctg in \$(seq 1 22) X; do
        loci=\$(ls *_loci_*_chr\${ctg}.txt 2>/dev/null | head -n1 || true)
        if [ -z "\$loci" ]; then
            echo "ERROR: no loci file matching *_loci_*_chr\${ctg}.txt among the staged ASCAT loci files" >&2
            exit 1
        fi
        ln -s ../../"\$loci" clairsto_cna_resources/loci_files/G1000_loci_hg38_chr\${ctg}.txt

        allele=\$(ls *_alleles_*_chr\${ctg}.txt 2>/dev/null | head -n1 || true)
        if [ -z "\$allele" ]; then
            echo "ERROR: no allele file matching *_alleles_*_chr\${ctg}.txt among the staged ASCAT allele files" >&2
            exit 1
        fi
        ln -s ../../"\$allele" clairsto_cna_resources/allele_files/G1000_alleles_hg38_chr\${ctg}.txt
    done

    # Verdict's correct_logr.py prepends "chr" to the contig column itself, so the GC table must
    # carry bare contig names. The ASCAT distribution used for the pipeline's own ASCAT step
    # carries them chr-prefixed; strip the prefix here rather than patching the container.
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{print; next} {sub(/^chr/, "", \$2); print}' \\
        ${gc_file} > clairsto_cna_resources/GC_G1000_hg38.txt

    # The replication timing table is generated against these same loci with bare contig names
    # (assets/clairsto_cna_chm13/make_chm13_rt.py), so it is used as-is.
    ln -s ../${rt_file} clairsto_cna_resources/RT_G1000_hg38.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/^.*Awk //; s/,.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clairsto_cna_resources/loci_files clairsto_cna_resources/allele_files
    for ctg in \$(seq 1 22) X; do
        touch clairsto_cna_resources/loci_files/G1000_loci_hg38_chr\${ctg}.txt
        touch clairsto_cna_resources/allele_files/G1000_alleles_hg38_chr\${ctg}.txt
    done
    touch clairsto_cna_resources/GC_G1000_hg38.txt
    touch clairsto_cna_resources/RT_G1000_hg38.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 5.1.0
    END_VERSIONS
    """
}
