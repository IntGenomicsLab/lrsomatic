process CLAIRSTO_VERDICT_TAG {
    tag "$meta.id"
    label 'process_low'

    // ClairS-TO's Verdict tagging step, run on ASCAT's purity and segments; same fork image as CLAIRSTO
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13-c0687e8-cpu'
        : 'ghcr.io/ljwharbers/clairs-to:0.5.1-verdict-chm13-c0687e8-cpu'}"

    input:
    // Untagged ClairS-TO VCFs (--disable_verdict) and ASCAT's tables; the VCFs keep their names, so they are staged aside
    tuple val(meta), path(snv_vcf, stageAs: 'untagged/snv.vcf.gz'), path(indel_vcf, stageAs: 'untagged/indel.vcf.gz'), path(purityploidy), path(segments)

    output:
    tuple val(meta), path("snv.vcf.gz"),                 emit: snv_vcf
    tuple val(meta), path("snv.vcf.gz.tbi"),             emit: snv_tbi
    tuple val(meta), path("indel.vcf.gz"),               emit: indel_vcf
    tuple val(meta), path("indel.vcf.gz.tbi"),           emit: indel_tbi
    // What the tags were computed from; absent when ASCAT found no solution
    tuple val(meta), path("*_Tumor_Purity_Ploidy.txt"),  emit: purity_ploidy, optional: true
    tuple val(meta), path("*_Tumor_CNA.txt"),            emit: cna
    tuple val("${task.process}"), val('clairsto'), eval("run_clairs_to --version |& sed '1!d ; s/run_clairs_to //'"), topic: versions, emit: versions_clairsto

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The image's tools live in its conda environment: python3 with scipy, bgzip, tabix
    export PATH=/opt/micromamba/envs/clairs-to/bin:\$PATH

    ascat_to_verdict.py \\
        --purityploidy $purityploidy \\
        --segments $segments \\
        --vcf untagged/snv.vcf.gz \\
        --sample ${prefix} \\
        --purity_out ${prefix}_Tumor_Purity_Ploidy.txt \\
        --cna_out ${prefix}_Tumor_CNA.txt

    for kind in snv indel; do
        if [ -s ${prefix}_Tumor_Purity_Ploidy.txt ]; then
            python3 /opt/bin/src/verdict/tag_germline_variant.py \\
                --input_vcf_fn untagged/\${kind}.vcf.gz \\
                --output_fn \${kind}.vcf \\
                --tumor_purity_ploidy_output_file ${prefix}_Tumor_Purity_Ploidy.txt \\
                --tumor_cna_output_file ${prefix}_Tumor_CNA.txt
        fi
        # Verdict's bgzip/tabix run unchecked and bgzip removes its input on success, so a leftover plain VCF means compression failed
        if [ -e \${kind}.vcf ]; then
            echo "ERROR: Verdict tagged \${kind} but compressing \${kind}.vcf failed" >&2
            exit 1
        fi

        # No tags (purity above 0.6, or no ASCAT solution): pass the calls through, as ClairS-TO does
        if [ ! -e \${kind}.vcf.gz ]; then
            cp untagged/\${kind}.vcf.gz \${kind}.vcf.gz
            tabix -f -p vcf \${kind}.vcf.gz
        fi
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > snv.vcf.gz
    touch snv.vcf.gz.tbi
    echo "" | gzip > indel.vcf.gz
    touch indel.vcf.gz.tbi
    printf 'Sample\\tPurity\\tPloidy\\tGoodnessOfFit\\n${prefix}\\t0.5\\t2.0\\tNA\\n' > ${prefix}_Tumor_Purity_Ploidy.txt
    printf 'Sample\\tChromosome\\tStartPosition\\tEndPosition\\tnMajor\\tnMinor\\n' > ${prefix}_Tumor_CNA.txt
    """
}
