process PROBE {
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: done

    script:
    """
    if [ ${task.attempt} -eq 1 ] || [ "${meta.always}" = true ]; then exit ${meta.exit}; fi
    echo "${task.process.tokenize(':').last()} ${task.attempt}" > ${meta.id}.txt
    """
}
