process HMF_LIFTOVER {
    tag "${chain.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0':
        'biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    // The GRCh38 HMF resources in LRSOMATICREPORT's slot order, `[]` for an absent optional one
    tuple path(driver_panel, stageAs: 'in/driver_panel/*'), path(ensembl_data, stageAs: 'in/ensembl/*'), path(known_fusions, stageAs: 'in/known_fusions/*'), path(hotspots, stageAs: 'in/hotspots/*'), path(actionable_bed, stageAs: 'in/actionable_bed/*')
    path(chain)                 // GRCh38 -> target chain (.chain / .chain.gz)
    tuple path(fasta), path(fai) // target FASTA, for the hotspot REF check and ##contig lines

    output:
    // The lifted resources laid out like an oncoanalyser bundle (common/, dna/sv/, dna/variants/), so the published
    // copy also works as --report_hmf_resource_dir; hmfLiftedFiles() maps it back to the report's slots
    path("hmf_chm13"), emit: bundle
    path("hmf_chm13/hmf_liftover_report.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opt = { flag, f -> f ? "${flag} ${f}" : '' }
    """
    hmf_liftover.py \\
        --chain ${chain} \\
        --fasta ${fasta} \\
        --outdir hmf_chm13 \\
        ${opt.call('--driver-panel', driver_panel)} \\
        ${opt.call('--ensembl-data-dir', ensembl_data)} \\
        ${opt.call('--known-fusions', known_fusions)} \\
        ${opt.call('--hotspots', hotspots)} \\
        ${opt.call('--actionable-bed', actionable_bed)} \\
        ${args}
    """

    stub:
    def touch = { dir, f -> f ? "touch hmf_chm13/${dir}/${f.name}${dir == 'dna/variants' && f.name.endsWith('.vcf') ? '.gz' : ''}" : '' }
    """
    mkdir -p hmf_chm13/common/ensembl_data hmf_chm13/dna/sv hmf_chm13/dna/variants
    echo t2t > hmf_chm13/common/ensembl_data/REFERENCE
    touch hmf_chm13/common/ensembl_data/ensembl_gene_data.csv hmf_chm13/common/ensembl_data/ensembl_trans_exon_data.csv
    ${touch.call('common', driver_panel)}
    ${touch.call('dna/sv', known_fusions)}
    ${touch.call('dna/variants', hotspots)}
    ${touch.call('dna/variants', actionable_bed)}
    printf 'resource\\titem\\tcount\\n' > hmf_chm13/hmf_liftover_report.tsv
    """
}
