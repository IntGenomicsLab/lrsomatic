process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    // No conda: the image ships render_report.R itself, not just its dependencies (guard in `script:`).
    // TODO: switch to bioconda `lrsomatic-report` once the recipe in ljwharbers/lrsomatic_report is merged. Version bump = these two tags.
    // TODO: pin to 2.0.0 once released
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/lrsomatic-report-sif:2.0.0rc1'
        : 'ghcr.io/ljwharbers/lrsomatic-report:2.0.0rc1'}"

    input:
    // Every path input is optional (`[]` when skipped); tumor/normal QC stage apart because a pair shares meta.id.
    // extra_relpaths/extra_files are parallel lists: extra_files[i] is linked at sample_dir/<extra_relpaths[i]>. Each file
    // stages in its own numbered dir (extra/1/, extra/2/, ...), so files sharing a base name cannot collide.
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), val(extra_relpaths), path(extra_files, stageAs: 'extra/?/*')
    // Builtin gene panel TSVs, owned by the pipeline; reach the tool as --gene-lists-dir
    path(gene_lists, stageAs: 'gene_lists')
    // User-supplied gene panel TSVs (`[]` for builtins); the `--gene-panel` args are built in conf/modules.config
    path(gene_panels, stageAs: 'gene_panels/*')

    output:
    tuple val(meta), path("*_report.html")          , emit: report
    // Written whenever a somatic VEP file is present
    tuple val(meta), path("*_small_variants.tsv.gz"), emit: small_variants, optional: true
    tuple val(meta), path("staged.txt")             , emit: staged        , optional: true // stub only: the sample_dir layout, for the tests
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
    def link_sample_dir = reportSampleDirCommands(vep_somatic, sv_vep, severus_vcf, somatic_vcf, ascat_files, extra_relpaths, extra_files)

    """
    # Quarto/Deno caches must live in the task dir: an inherited HOME/TMPDIR/XDG_CACHE_HOME is read-only inside the container
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$TMPDIR"

    ${link_sample_dir}

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
    // Builds the same sample_dir as the script, so a -stub run shows what the tool would discover
    def link_sample_dir = reportSampleDirCommands(vep_somatic, sv_vep, severus_vcf, somatic_vcf, ascat_files, extra_relpaths, extra_files)
    def touch_small_variants = vep_somatic ? "echo '' | gzip > ${prefix}_small_variants.tsv.gz" : ''
    """
    ${link_sample_dir}

    find -L sample_dir -type f | LC_ALL=C sort | while read -r f; do
        printf '%s\\t%s\\n' "\$f" "\$(md5sum < "\$f" | cut -d' ' -f1)"
    done > staged.txt
    touch ${prefix}_report.html
    ${touch_small_variants}
    """
}

// Shell commands that build sample_dir/, mirroring the published <outdir>/<id>/ layout the tool discovers recursively.
// Links go file by file: R's list.files(recursive = TRUE) does not descend into symlinked dirs.
def reportSampleDirCommands(vep_somatic, sv_vep, severus_vcf, somatic_vcf, ascat_files, extra_relpaths, extra_files) {
    // Discovery matches on base name, so suffix-distinct files can be linked flat
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs
        ? "for f in ${flat_inputs.collect { f -> "\"${f}\"" }.join(' ')}; do ln -s \"\$PWD/\$f\" \"sample_dir/\$f\"; done"
        : ''

    // The VAF/depth/phasing source is looked up at a literal path
    def link_somatic = somatic_vcf
        ? "mkdir -p sample_dir/variants/phased && ln -s \"\$PWD/${somatic_vcf}\" sample_dir/variants/phased/somatic_smallvariants.vcf.gz"
        : ''

    // Every other output lands at its relpath under sample_dir; the two lists pair up by position
    def extra_staged = [extra_files].flatten().findAll { f -> f }
    def extra_rel    = [extra_relpaths].flatten().findAll { r -> r }.collect { r -> r.toString() }
    if (extra_staged.size() != extra_rel.size()) {
        error "LRSOMATICREPORT: ${extra_rel.size()} extra_relpaths for ${extra_staged.size()} extra_files"
    }
    if (extra_rel.any { r -> r.startsWith('/') || r.tokenize('/').contains('..') } || extra_rel.unique(false).size() != extra_rel.size()) {
        error "LRSOMATICREPORT: extra_relpaths must be unique and relative: ${extra_rel}"
    }
    def link_extra = [extra_rel, extra_staged].transpose().collect { rel, f ->
        def dest = "sample_dir/${rel}"
        "mkdir -p \"${dest.substring(0, dest.lastIndexOf('/'))}\" && ln -s \"\$PWD/${f}\" \"${dest}\""
    }.join('\n    ')

    return """\
    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi
    ${link_extra}"""
}
