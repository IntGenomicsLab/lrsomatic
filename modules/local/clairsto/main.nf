process CLAIRSTO {
    tag "$meta.id"
    label 'process_very_high'

    // Patched build of ClairS-TO v0.5.1. The patch is a three-line fix to
    // src/verdict/run_ascat.py, which raises an IndexError when a logR segment has no
    // heterozygous SNP within the +/-10000-probe window (hit on the single-copy part of a male
    // chrX under CHM13); it falls back to the nearest heterozygous probe. The file is byte
    // identical from v0.4.2 to v0.5.1 upstream, so this still applies to HEAD.
    // Return to docker.io/hkubal/clairs-to once the fix is released upstream.
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/ljwharbers/clairs-to-sif:0.5.1-verdictchm13-REPLACEME':
        'ghcr.io/ljwharbers/clairs-to:0.5.1-verdictchm13-REPLACEME' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model), path(pon_vcfs), val(pon_flags)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)
    path(cna_resources)

    output:
    tuple val(meta), path("indel.vcf.gz"),      emit: indel_vcf
    tuple val(meta), path("indel.vcf.gz.tbi"),  emit: indel_tbi
    tuple val(meta), path("snv.vcf.gz"),        emit: snv_vcf
    tuple val(meta), path("snv.vcf.gz.tbi"),    emit: snv_tbi
    tuple val(meta), path("*_verdict.log"),               optional: true, emit: verdict_log
    tuple val(meta), path("*_verdict_purity_ploidy.txt"), optional: true, emit: verdict_purity
    tuple val("${task.process}"), val('clairsto'), eval("run_clairs_to  --version |& sed '1!d ; s/run_clairs_to //'"), topic: versions, emit: versions_clairsto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def conda_prefix = (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''
    def pon_string   = pon_vcfs.join(',')
    def flags_string = pon_flags.join(',')
    // Verdict resource set matching the reference build; empty means the GRCh38 set inside the image
    def cna = cna_resources ? "--cna_resource_dir ${cna_resources}" : ''

    // Verdict never fails loudly: every failure path inside it exits 0, and on a non-GRCh38
    // reference it silently tags nothing (its contig check only looks at chr1..chrX naming).
    // Check afterwards that it actually ran, at the severity set by --clairsto_verdict_check.
    def verdict_check = task.ext.args2 ?: 'off'
    def check_verdict = verdict_check != 'off' && !args.contains('--disable_verdict')
    def level = verdict_check == 'warn' ? 'WARNING' : 'ERROR'
    def fail  = verdict_check == 'warn' ? '' : '\n        exit 1'
    def verdict_guard = !check_verdict ? '' : """
    # Verdict sanity check (--clairsto_verdict_check ${verdict_check})
    run_log=\$(ls run_clairs_to*.log 2>/dev/null | head -n1 || true)
    cgt_log=\$(ls logs*/6_CGT.log 2>/dev/null | head -n1 || true)
    verdict_reason=''
    if [ -n "\$run_log" ] && grep -q 'DISABLE APPLYING VERDICT: True' "\$run_log"; then
        verdict_reason='ClairS-TO disabled Verdict during its pre-flight checks'
    elif [ -z "\$cgt_log" ]; then
        verdict_reason='ClairS-TO produced no Verdict log (logs*/6_CGT.log)'
    elif grep -q 'ERROR in STEP' "\$cgt_log"; then
        verdict_reason='a Verdict sub-step failed (see the published Verdict log)'
    elif grep -q 'Verdict can not obtain final results' "\$cgt_log"; then
        verdict_reason='Verdict could not estimate purity/ploidy or copy-number segments'
    elif [ "\$(zcat snv.vcf.gz | grep -c 'Verdict_' || true)" -eq 0 ]; then
        verdict_reason='Verdict tagged 0 variants'
    fi
    if [ -n "\$verdict_reason" ]; then
        echo "[lrsomatic] ${level}: ClairS-TO germline tagging (Verdict) did not run for ${prefix}: \$verdict_reason." >&2
        echo "[lrsomatic] Verdict is what removes high-VAF germline variants from a tumor-only call set; without it they stay PASS." >&2
        echo "[lrsomatic] Pass --clairsto_cna_resources <dir> with a resource set built for this reference," >&2
        echo "[lrsomatic] or accept the loss with --clairsto_disable_verdict, or downgrade this to a warning with --clairsto_verdict_check warn." >&2${fail}
    fi
"""

    """
    /opt/bin/run_clairs_to \\
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
        $cna \\
        $args

    # From v0.4.4 ClairS-TO rewrites a *default* output prefix whenever --sample_name is not
    # "SAMPLE" (snv -> snv_<sample>), and passing `--snv_output_prefix snv` does not opt out
    # because the guard compares against the literal "snv". Ask for prefixes it leaves alone,
    # then restore the canonical names the rest of the pipeline expects.
    mv snv_out.vcf.gz snv.vcf.gz
    mv snv_out.vcf.gz.tbi snv.vcf.gz.tbi
    mv indel_out.vcf.gz indel.vcf.gz
    mv indel_out.vcf.gz.tbi indel.vcf.gz.tbi
${verdict_guard}
    # Keep the Verdict log and the purity/ploidy estimate as QC. Both live under work directories
    # whose names are sample-scoped from v0.4.4 on (tmp_<sample>/, logs_<sample>/), hence the globs.
    # The purity file is named from --tumor_sample_name (default "tumor"), so rename it per sample.
    qc_log=\$(ls logs*/6_CGT.log 2>/dev/null | head -n1 || true)
    if [ -n "\$qc_log" ]; then cp "\$qc_log" ${prefix}_verdict.log; fi
    qc_pp=\$(ls tmp*/cna_output/*_Tumor_Purity_Ploidy.txt 2>/dev/null | head -n1 || true)
    if [ -n "\$qc_pp" ]; then cp "\$qc_pp" ${prefix}_verdict_purity_ploidy.txt; fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p output
    echo "" | gzip > snv.vcf.gz
    touch snv.vcf.gz.tbi
    echo "" | gzip > indel.vcf.gz
    touch indel.vcf.gz.tbi
    touch ${prefix}_verdict.log
    touch ${prefix}_verdict_purity_ploidy.txt
    """
}
