process CLAIRSTO_CNA_RESOURCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    // The same ASCAT loci/allele/GC set the ASCAT subworkflow uses, already unzipped
    tuple val(meta), path(loci, stageAs: 'loci/*'), path(alleles, stageAs: 'alleles/*'), path(gc, stageAs: 'gc/*')

    output:
    tuple val(meta), path("cna_resources"), emit: cna_resources
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p cna_resources/loci_files cna_resources/allele_files

    # Recursive: a set unzipped from an archive arrives as one directory holding the
    # per-contig files, while one given as a plain path arrives as the files themselves.
    # Real files either way -- only this directory is staged into CLAIRSTO, so links out of
    # it would dangle.
    find -L loci -type f -exec cp -L -t cna_resources/loci_files/ {} +
    find -L alleles -type f -exec cp -L -t cna_resources/allele_files/ {} +

    # ClairS-TO takes exactly one GC_*.txt and errors on a second candidate. Fail here instead,
    # where the message can name the files.
    n_gc=\$(find -L gc -type f | wc -l)
    if [ "\$n_gc" -ne 1 ]; then
        echo "ERROR: expected exactly one GC content file, found \$n_gc:" >&2
        find -L gc -type f >&2
        exit 1
    fi

    gc_src=\$(find -L gc -type f)
    case "\$(basename "\$gc_src")" in
        GC_*.txt) cp -L "\$gc_src" "cna_resources/\$(basename "\$gc_src")" ;;
        *)        cp -L "\$gc_src" "cna_resources/GC_${prefix}.txt" ;;
    esac

    # No RT_*.txt, deliberately: Verdict then corrects LogR for GC content only, which is what
    # ClairS-TO recommends for CHM13.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cp --version | sed '1!d; s/cp (GNU coreutils) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p cna_resources/loci_files cna_resources/allele_files
    for contig in \$(seq 1 22) X; do
        touch "cna_resources/loci_files/G1000_loci_${prefix}_chr\${contig}.txt"
        touch "cna_resources/allele_files/G1000_alleles_${prefix}_chr\${contig}.txt"
    done
    touch "cna_resources/GC_${prefix}.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: 9.5
    END_VERSIONS
    """
}
