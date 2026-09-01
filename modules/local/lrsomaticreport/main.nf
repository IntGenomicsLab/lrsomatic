process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Dependencies only (R, Quarto and the tool's R packages), built via the Wave
    // containers API from this module's environment.yml (frozen build). The tool
    // itself is vendored at assets/lrsomatic_report -- see VENDORED.md there.
    // Two separate Wave builds are needed: `wave --singularity` produces a
    // Singularity-native SIF artifact (the oras:// reference), while the default
    // build produces a genuine OCI image (the plain tag). Rebuild both whenever
    // environment.yml changes:
    //   wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await [--singularity]
    // NOTE (local patch, not upstream): dropped the `&& !task.ext.singularity_pull_docker_container`
    // guard -- `task` is not in scope inside a `container` directive closure on this cluster's
    // Nextflow build (25.09.0-beta), so the original line threw "No such variable: task" and
    // killed the whole run at the very last step. Nothing in this project's config ever sets
    // ext.singularity_pull_docker_container for this process, so the guard was always a no-op
    // for us anyway -- safe to drop rather than work around.
    container "${workflow.containerEngine == 'singularity'
        ? 'oras://community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:dc62d809aa6fd497'
        : 'community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:c1049dbaf31bf178'}"

    input:
    // All per-sample report inputs are optional (path may be `[]` if the corresponding
    // upstream tool was skipped or produced no output for this sample); the report tool
    // renders a "not available" notice for any missing section.
    // qc_tumor_files/qc_normal_files are staged into distinct subdirectories:
    // mosdepth/samtools default to a `${meta.id}`-only prefix (see conf/modules.config),
    // so for a matched T/N pair (same meta.id) the tumor and normal QC files are
    // identically named -- staging both lists flat would collide.
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), path(wakhan_files, stageAs: 'wakhan/*')
    path(report_src) // lrsomatic_report source tree (bin/, R/, templates/, assets/)
    // Optional user-supplied gene panel TSV, staged so it is bound into the container;
    // `[]` when --report_gene_panel names a builtin panel (or is unset), which the tool
    // resolves from report_src/assets/gene_lists instead. Either way the `--gene-panel`
    // argument itself is built in conf/modules.config, which passes the *base* name --
    // the staged name of this file when it is one.
    path(gene_panel)

    output:
    tuple val(meta), path("*_report.html"), emit: report
    // No CLI version flag is provided by the tool; keep in sync with the vendored
    // release recorded in assets/lrsomatic_report/VENDORED.md
    tuple val("${task.process}"), val('lrsomatic_report'), val("1.2.1"), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'

    // Discovery (R/locate_outputs.R, R/sections/sv.R) is recursive under sample_dir and
    // matches on the *base name*, so anything identified by a distinctive filename suffix
    // can be linked flat: the VEP somatic VCF (*_SOMATIC_VEP.vcf.gz), the VEP SV VCF
    // (*_SV_VEP.vcf.gz), the Severus SV VCF (severus_somatic.vcf.gz) and every ASCAT file
    // (*.segments_raw.txt, *.purityploidy.txt, the diagnostic PNGs).
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs ? """
    for f in ${flat_inputs.join(' ')}; do ln -s "\$PWD/\$f" "sample_dir/\$f"; done
    """ : ''

    // The VAF/depth/phasing source is the exception: locate_outputs() looks for it at the
    // literal path variants/phased/somatic_smallvariants.vcf.gz before falling back to a
    // caller-specific directory, so this one file needs its canonical name and location.
    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/variants/phased
    ln -s "\$PWD/${somatic_vcf}" sample_dir/variants/phased/somatic_smallvariants.vcf.gz
    """ : ''

    """
    # Quarto/Deno write a cache dir under \$HOME and a session dir under \$TMPDIR;
    # point both at the task work dir, which is always writable and always bound into
    # the container. Relying on the container's own \$HOME and /tmp fails on clusters
    # that mount them read-only ("Read-only file system (os error 30): tmpdir").
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp
    mkdir -p "\$TMPDIR"

    # The Wave/conda-built container doesn't auto-source conda's activation hooks
    # (e.g. quarto needs QUARTO_SHARE_PATH); source them if present. Some hooks
    # (e.g. gcc_linux-64) reference \$CONDA_PREFIX under `set -u`, so export it first.
    # Under `-profile conda` CONDA_PREFIX already points at the task's own env, so
    # only fall back to the container's /opt/conda when it is unset.
    export CONDA_PREFIX="\${CONDA_PREFIX:-/opt/conda}"
    for f in "\$CONDA_PREFIX"/etc/conda/activate.d/*.sh; do
        [ -f "\$f" ] && source "\$f"
    done

    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}

    # QC is split tumor/normal by path: locate_outputs() takes the first suffix match
    # outside any /normal/ component as tumor, and the first one inside it as normal.
    # Link file by file rather than symlinking the staging directory itself -- R's
    # list.files(recursive = TRUE) does not descend into symlinked directories.
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi

    # Wakhan is addressed by fixed path, not by suffix: sample_dir/wakhan must hold
    # solutions_ranks.tsv, *heatmap_ploidy_purity.html and the per-solution
    # solution_<rank>/ directories (R/parse_ascat.R:locate_wakhan_cn_plots).
    if [ -d wakhan ]; then
        mkdir -p sample_dir/wakhan
        for f in wakhan/*; do ln -s "\$PWD/\$f" "sample_dir/wakhan/\$(basename "\$f")"; done
    fi

    Rscript ${report_src}/bin/render_report.R \\
        --sample-dir sample_dir \\
        --sample-id ${prefix} \\
        --sex ${sex} \\
        --reference auto \\
        --output ${prefix}_report.html \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
