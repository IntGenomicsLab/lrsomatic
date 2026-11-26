process HMF_BUNDLE {
    tag "${url.toString().tokenize('/').last()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0':
        'biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    // The genome's hosted HMF report bundle (conf/igenomes.config hmf_bundle_url / hmf_bundle_md5); `val` so the
    // task hash, and with it -resume, depends on the URL and MD5 only
    tuple val(url), val(md5)

    output:
    // Unpacked without its top-level directory: common/, dna/sv/, dna/variants/ (hmfHostedFiles() maps the slots)
    path("hmf_bundle"), emit: bundle
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fetch_hmf_bundle.py \\
        --url '${url}' \\
        --md5 '${md5 ?: ''}' \\
        --outdir hmf_bundle
    """

    stub:
    """
    mkdir -p hmf_bundle/common/ensembl_data hmf_bundle/dna/sv hmf_bundle/dna/variants
    touch hmf_bundle/common/DriverGenePanel.38.tsv
    touch hmf_bundle/common/ensembl_data/ensembl_gene_data.csv hmf_bundle/common/ensembl_data/ensembl_trans_exon_data.csv
    touch hmf_bundle/dna/sv/known_fusion_data.38.csv
    touch hmf_bundle/dna/variants/KnownHotspots.somatic.38.vcf.gz hmf_bundle/dna/variants/ActionableCodingPanel.38.bed.gz
    """
}
