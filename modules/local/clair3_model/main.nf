process CLAIR3_MODEL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data'
        : 'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e'}"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${meta.id}", type: 'dir'), emit: model
    tuple val("${task.process}"), val('wget'), eval('wget --version | head -1 | cut -d " " -f 3'), emit: versions_wget, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Clair3 v2 models are directories holding two PyTorch checkpoints, served unpacked. Clair3 unpickles
    // them, so the certificate is checked: wget's default CA path is empty in the image, but conda ships a bundle
    """
    mkdir ${meta.id}
    for checkpoint in pileup.pt full_alignment.pt; do
        wget ${args} --ca-certificate=\${CONDA_PREFIX:-/opt/conda}/ssl/cert.pem -O ${meta.id}/\$checkpoint ${url}/\$checkpoint
    done
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/pileup.pt ${meta.id}/full_alignment.pt
    """
}
