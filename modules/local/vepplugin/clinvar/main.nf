process VEPPLUGIN_CLINVAR {
    tag "${vcf_url.toString().tokenize('/').last()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data'
        : 'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e'}"

    input:
    tuple val(vcf_url), val(tbi_url), val(md5), val(tbi_md5)

    output:
    path "${vcf_name}{,.tbi}", emit: files
    // versions.yml rather than an eval() topic: eval outputs are numbered pipeline-wide, so adding
    // one shifts the cache key of every other task that has one and breaks -resume of existing runs
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    vcf_name = vcf_url.toString().tokenize('/').last()
    def sums = [ md5 ? "${md5}  ${vcf_name}" : null, tbi_md5 ? "${tbi_md5}  ${vcf_name}.tbi" : null ].findAll()
    def check = sums ? "printf '%s\\n' ${sums.collect { line -> "'${line}'" }.join(' ')} | md5sum -c -" : ''
    """
    # Retry a transient 503 from NCBI: a linear backoff of 1 s, 2 s, ... up to 10 s, about 45 s over 10 tries
    wget \\
        --no-verbose \\
        --tries=10 \\
        --waitretry=10 \\
        --retry-on-http-error=429,500,502,503,504 \\
        ${args} \\
        -O ${vcf_name} \\
        ${vcf_url}

    # Saved next to the VCF under the name VEP looks for, whatever the host calls it
    wget \\
        --no-verbose \\
        --tries=10 \\
        --waitretry=10 \\
        --retry-on-http-error=429,500,502,503,504 \\
        ${args} \\
        -O ${vcf_name}.tbi \\
        ${tbi_url}

    # Pinned checksums keep the release fixed: a host that re-publishes under the same name fails here
    ${check}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    vcf_name = vcf_url.toString().tokenize('/').last()
    """
    echo "" | gzip > ${vcf_name}
    touch ${vcf_name}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
