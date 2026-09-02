process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Dependencies only -- the tool itself is vendored at assets/lrsomatic_report.
    // Singularity needs a Wave-built SIF (oras://), Docker a plain OCI image; rebuild
    // both whenever environment.yml changes:
    //   wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await [--singularity]
    container "${workflow.containerEngine == 'singularity'
        ? 'oras://community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:dc62d809aa6fd497'
        : 'community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:c1049dbaf31bf178'}"

    input:
    // Every path input is optional (`[]` when the upstream tool was skipped); the report
    // renders a "not available" notice for the missing section. Tumor and normal QC stage
    // into separate subdirectories because a matched pair shares meta.id, so their QC
    // files are identically named and would collide if staged flat.
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), path(wakhan_files, stageAs: 'wakhan/*')
    path(report_src) // lrsomatic_report source tree (bin/, R/, templates/, assets/)
    // Optional user-supplied gene panel TSVs, staged so they are bound into the container;
    // `[]` for builtin panels, which the tool resolves from report_src itself. The matching
    // `--gene-panel gene_panels/<base>` arguments are built in conf/modules.config.
    path(gene_panels, stageAs: 'gene_panels/*')

    output:
    tuple val(meta), path("*_report.html"), emit: report
    // WARN: Manually update to match the vendored release in assets/lrsomatic_report/VENDORED.md
    tuple val("${task.process}"), val('lrsomatic_report'), val('1.3.0'), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'

    // Discovery is recursive under sample_dir and matches on the base name, so anything
    // with a distinctive filename suffix can be linked flat
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs ? """
    for f in ${flat_inputs.collect { f -> "\"${f}\"" }.join(' ')}; do ln -s "\$PWD/\$f" "sample_dir/\$f"; done
    """ : ''

    // The exception: the VAF/depth/phasing source is looked up at a literal path, so it
    // needs its canonical name and location
    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/variants/phased
    ln -s "\$PWD/${somatic_vcf}" sample_dir/variants/phased/somatic_smallvariants.vcf.gz
    """ : ''

    """
    # Quarto/Deno write a cache dir under \$HOME and a session dir under \$TMPDIR; point
    # both at the task work dir, as clusters may mount the container's \$HOME and /tmp
    # read-only
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp
    mkdir -p "\$TMPDIR"

    # The Wave/conda-built container doesn't auto-source conda's activation hooks (quarto
    # needs QUARTO_SHARE_PATH). Some hooks reference \$CONDA_PREFIX under `set -u`, so export
    # it first -- falling back to the container's /opt/conda only when unset, as
    # `-profile conda` already points it at the task's own env
    export CONDA_PREFIX="\${CONDA_PREFIX:-/opt/conda}"
    for f in "\$CONDA_PREFIX"/etc/conda/activate.d/*.sh; do
        [ -f "\$f" ] && source "\$f"
    done

    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}

    # QC is split tumor/normal by path. Link file by file rather than symlinking the
    # staging directory -- R's list.files(recursive = TRUE) does not descend into symlinks
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi

    # Wakhan is addressed by fixed path, not by suffix: sample_dir/wakhan must hold
    # solutions_ranks.tsv, the heatmap HTML and the per-solution solution_<rank>/ directories
    if [ -d wakhan ]; then
        mkdir -p sample_dir/wakhan
        for f in wakhan/*; do ln -s "\$PWD/\$f" "sample_dir/wakhan/\$(basename "\$f")"; done
    fi

    Rscript "${report_src}/bin/render_report.R" \\
        --sample-dir sample_dir \\
        --sample-id "${prefix}" \\
        --sex "${sex}" \\
        --reference auto \\
        --output "${prefix}_report.html" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
