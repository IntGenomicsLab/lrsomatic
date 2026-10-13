process AADATAREPO_DOWNLOAD {
    tag "${url.toString().tokenize('/').last()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data'
        : 'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e'}"

    input:
    tuple val(meta), val(url), val(md5)

    output:
    tuple val(meta), path("${archive_name}"), emit: archive
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    archive_name = url.toString().tokenize('/').last()
    def check = md5 ? "echo '${md5}  ${archive_name}' | md5sum -c -" : ''
    """
    wget \\
        --no-verbose \\
        --tries=10 \\
        --waitretry=10 \\
        ${args} \\
        -O ${archive_name} \\
        ${url}

    # The tarball is updated in place, so a pinned MD5 fails the task on a re-published repo
    ${check}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    archive_name = url.toString().tokenize('/').last()
    """
    echo "" | gzip > ${archive_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
