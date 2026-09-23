process CLAIRSTO {
    tag "$meta.id"
    label 'process_very_high'

    // Fork of ClairS-TO 0.5.1 that resolves Verdict's CNA resources from --cna_resource_dir
    // instead of hardcoded GRCh38 names. No conda build; revert once upstream carries it.
    // Hosted on Docker Hub, not ghcr: ghcr redirects blob downloads to an Azure URL that expires at the
    // next 5-minute mark and resets a stream still open then, and Apptainer resumes neither an oras://
    // nor a docker:// download, so this 3.3 GB SIF failed on any link slower than ~10 MB/s. Docker Hub's
    // CloudFront URL is valid for 50 minutes and is only checked when the request starts.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13-c0687e8-flat'
        : 'docker.io/ljwharbers/clairs-to:0.5.1-verdict-chm13-c0687e8-flat'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model), path(pon_vcfs), val(pon_flags)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)
    // Verdict's ASCAT set for this assembly, from CLAIRSTO_CNA_RESOURCES; [] uses the image's own
    tuple val(meta4), path(cna_resources)

    output:
    tuple val(meta), path("indel.vcf.gz"),      emit: indel_vcf
    tuple val(meta), path("indel.vcf.gz.tbi"),  emit: indel_tbi
    tuple val(meta), path("snv.vcf.gz"),        emit: snv_vcf
    tuple val(meta), path("snv.vcf.gz.tbi"),    emit: snv_tbi
    // What Verdict tagged from, when it ran inside ClairS-TO (--skip_ascat); absent otherwise
    tuple val(meta), path("*_Tumor_Purity_Ploidy.txt"), emit: purity_ploidy, optional: true
    tuple val(meta), path("*_Tumor_CNA.txt"),           emit: cna,           optional: true
    tuple val("${task.process}"), val('clairsto'), eval("run_clairs_to  --version |& sed '1!d ; s/run_clairs_to //'"), topic: versions, emit: versions_clairsto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The launcher cannot locate its own conda environment inside a read-only image
    def conda_prefix = (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''
    // Omitted for GRCh38; a set that cannot belong to --ref_fn disables Verdict with a warning
    def cna_resource_dir = cna_resources ? "--cna_resource_dir ${cna_resources}" : ''
    def pon_string   = pon_vcfs.join(',')
    def flags_string = pon_flags.join(',')

    """
    /opt/bin/run_clairs_to \
        --tumor_bam_fn $tumor_bam \\
        --ref_fn $reference \\
        --platform $model \\
        --threads $task.cpus \\
        --output_dir . \\
        --sample_name ${prefix} \\
        --snv_output_prefix snv_out \\
        --indel_output_prefix indel_out \\
        --panel_of_normals ${pon_string} \\
        --panel_of_normals_require_allele_matching ${flags_string} \\
        $conda_prefix \\
        $cna_resource_dir \\
        $args

    # Explicit prefixes are left alone by 0.5.1's --sample_name renaming, so these names are fixed
    mv snv_out.vcf.gz snv.vcf.gz
    mv snv_out.vcf.gz.tbi snv.vcf.gz.tbi
    mv indel_out.vcf.gz indel.vcf.gz
    mv indel_out.vcf.gz.tbi indel.vcf.gz.tbi

    # Lift Verdict's own purity/CN tables out of ClairS-TO's work dir under the names CLAIRSTO_VERDICT_TAG publishes
    for table in Purity_Ploidy CNA; do
        src=\$(find . -path "*/cna_output/*_Tumor_\${table}.txt" -print -quit)
        if [ -n "\$src" ]; then
            cp -- "\$src" "${prefix}_Tumor_\${table}.txt"
        fi
    done
    """

    stub:
    """
    mkdir -p output
    echo "" | gzip > snv.vcf.gz
    touch snv.vcf.gz.tbi
    echo "" | gzip > indel.vcf.gz
    touch indel.vcf.gz.tbi
    """
}
