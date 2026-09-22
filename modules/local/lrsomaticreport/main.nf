process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    // No conda: the image ships render_report.R itself, not just its dependencies (guard in `script:`).
    // TODO: switch to bioconda `lrsomatic-report` once the recipe in ljwharbers/lrsomatic_report is merged. Version bump = these two tags.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/lrsomatic-report-sif:1.6.0'
        : 'ghcr.io/ljwharbers/lrsomatic-report:1.6.0'}"

    input:
    // Every path input is optional (`[]` when skipped); tumor/normal QC stage apart because a pair shares meta.id
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), path(wakhan_files, stageAs: 'wakhan/*')
    // Builtin gene panel TSVs, owned by the pipeline; reach the tool as --gene-lists-dir
    path(gene_lists, stageAs: 'gene_lists')
    // User-supplied gene panel TSVs (`[]` for builtins); the `--gene-panel` args are built in conf/modules.config
    path(gene_panels, stageAs: 'gene_panels/*')

    output:
    tuple val(meta), path("*_report.html"), emit: report
    // `env -u R_HOME`: Apptainer forwards the host env, and R's R_HOME warning goes to stdout, polluting the version string
    tuple val("${task.process}"), val('lrsomatic_report'), eval('env -u R_HOME render_report.R --version'), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "LRSOMATICREPORT does not support Conda: the report tool ships only inside its container. Use Docker / Singularity / Apptainer, or --skip_report."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'

    // Discovery is recursive and matches on base name, so suffix-distinct files can be linked flat
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs ? """
    for f in ${flat_inputs.collect { f -> "\"${f}\"" }.join(' ')}; do ln -s "\$PWD/\$f" "sample_dir/\$f"; done
    """ : ''

    // The VAF/depth/phasing source is looked up at a literal path
    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/variants/phased
    ln -s "\$PWD/${somatic_vcf}" sample_dir/variants/phased/somatic_smallvariants.vcf.gz
    """ : ''

    """
    # Quarto/Deno caches must live in the task dir: an inherited HOME/TMPDIR/XDG_CACHE_HOME is read-only inside the container
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$TMPDIR"

    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}

    # Link file by file: R's list.files(recursive = TRUE) does not descend into symlinked dirs
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi

    # Wakhan is read from a fixed path: sample_dir/wakhan with solutions_ranks.tsv, the heatmap and solution_<rank>/
    if [ -d wakhan ]; then
        mkdir -p sample_dir/wakhan
        for f in wakhan/*; do ln -s "\$PWD/\$f" "sample_dir/wakhan/\$(basename "\$f")"; done
    fi

    render_report.R \\
        --sample-dir sample_dir \\
        --sample-id "${prefix}" \\
        --sex "${sex}" \\
        --reference auto \\
        --gene-lists-dir gene_lists \\
        --output "${prefix}_report.html" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
