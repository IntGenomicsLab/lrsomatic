process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built via the Wave containers API from this module's environment.yml (frozen
    // build). Two separate Wave builds were needed: a singularity.enabled=true
    // session only produces a Singularity-native SIF artifact (blob URL below),
    // while a docker.enabled=true session produces a genuine OCI image (plain tag).
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e0d4fabb2f79dcc0d3446f1bda84507eb52ac21ebea75fd29ee5b1b26c61ee34/data'
        : 'community.wave.seqera.io/library/r-base_quarto_r-data.table_r-dplyr_pruned:9d12b9297c3c4d38'}"

    input:
    // All per-sample report inputs are optional (path may be `[]` if the corresponding
    // upstream tool was skipped or produced no output for this sample); the report tool
    // renders a "not available" notice for any missing section.
    // qc_tumor_files/qc_normal_files are staged into distinct subdirectories:
    // mosdepth/samtools default to a `${meta.id}`-only prefix (see conf/modules.config),
    // so for a matched T/N pair (same meta.id) the tumor and normal QC files are
    // identically named -- staging both lists flat would collide.
    tuple val(meta), path(vep_somatic), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*')
    path(report_src) // staged lrsomatic_report repo (bin/, R/, templates/, assets/)

    output:
    tuple val(meta), path("*_report.html"), emit: report
    // No CLI version flag is provided by the tool; footer literal is "LRSomatic report v1.0" (templates/per_sample.qmd)
    tuple val("${task.process}"), val('lrsomatic_report'), val("1.0"), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'
    // matched (T/N) samples publish under variants/clairs/; tumor-only samples under variants/clairsto/
    // -- this only controls the report tool's run-mode detection/labelling, see locate_outputs.R
    def somatic_dir = meta.paired_data ? 'variants/clairs' : 'variants/clairsto'

    def link_vep = vep_somatic ? """
    mkdir -p sample_dir/vep/somatic
    ln -s "\$PWD/${vep_somatic}" "sample_dir/vep/somatic/${prefix}_SOMATIC_VEP.vcf.gz"
    """ : ''

    def link_severus = severus_vcf ? """
    mkdir -p sample_dir/variants/severus/somatic_SVs
    ln -s "\$PWD/${severus_vcf}" "sample_dir/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
    """ : ''

    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/${somatic_dir}
    ln -s "\$PWD/${somatic_vcf}" "sample_dir/${somatic_dir}/somatic.vcf.gz"
    """ : ''

    def ascat_file_list = ascat_files ? ascat_files.join(' ') : ''
    def link_ascat = ascat_files ? """
    mkdir -p sample_dir/ascat
    for f in ${ascat_file_list}; do ln -s "\$PWD/\$f" "sample_dir/ascat/\$f"; done
    """ : ''

    // $f includes the 'qc_tumor/' staging subdirectory (see stageAs above); the
    // destination link name uses just the basename.
    def qc_tumor_file_list = qc_tumor_files ? qc_tumor_files.join(' ') : ''
    def link_qc_tumor = qc_tumor_files ? """
    mkdir -p sample_dir/qc/tumor/mosdepth sample_dir/qc/tumor/cramino_aln sample_dir/qc/tumor/samtools
    for f in ${qc_tumor_file_list}; do
        fname=\$(basename "\$f")
        case "\$fname" in
            *.mosdepth.*.txt)    ln -s "\$PWD/\$f" "sample_dir/qc/tumor/mosdepth/\$fname" ;;
            *_cramino.txt)       ln -s "\$PWD/\$f" "sample_dir/qc/tumor/cramino_aln/\$fname" ;;
            *.flagstat|*.stats)  ln -s "\$PWD/\$f" "sample_dir/qc/tumor/samtools/\$fname" ;;
        esac
    done
    """ : ''

    def qc_normal_file_list = qc_normal_files ? qc_normal_files.join(' ') : ''
    def link_qc_normal = qc_normal_files ? """
    mkdir -p sample_dir/qc/normal/mosdepth sample_dir/qc/normal/cramino_aln sample_dir/qc/normal/samtools
    for f in ${qc_normal_file_list}; do
        fname=\$(basename "\$f")
        case "\$fname" in
            *.mosdepth.*.txt)    ln -s "\$PWD/\$f" "sample_dir/qc/normal/mosdepth/\$fname" ;;
            *_cramino.txt)       ln -s "\$PWD/\$f" "sample_dir/qc/normal/cramino_aln/\$fname" ;;
            *.flagstat|*.stats)  ln -s "\$PWD/\$f" "sample_dir/qc/normal/samtools/\$fname" ;;
        esac
    done
    """ : ''

    """
    # Quarto/Deno write a cache dir under \$HOME; point it at the task work dir
    # (always writable) rather than relying on the container's \$HOME being bound.
    export HOME=\$PWD

    # The Wave/conda-built container doesn't auto-source conda's activation hooks
    # (e.g. quarto needs QUARTO_SHARE_PATH); source them if present. Some hooks
    # (e.g. gcc_linux-64) reference \$CONDA_PREFIX under `set -u`, so export it first.
    export CONDA_PREFIX=/opt/conda
    for f in /opt/conda/etc/conda/activate.d/*.sh; do
        [ -f "\$f" ] && source "\$f"
    done

    mkdir -p sample_dir
    ${link_vep}
    ${link_severus}
    ${link_somatic}
    ${link_ascat}
    ${link_qc_tumor}
    ${link_qc_normal}

    # Quarto renders in-place next to the .qmd it's given (render_report.R's own
    # post-render step relies on this). report_src is a single fixed path shared
    # by every sample's task, so dereferencing it into a private, task-local copy
    # avoids concurrent per-sample renders colliding on the same physical directory.
    cp -rL "${report_src}" report_src_local

    Rscript report_src_local/bin/render_report.R \\
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
