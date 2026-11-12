// Test helpers, not part of the pipeline: build the inputs of the PHASING_HAPLOTYPING test from the CHM13 chr19 HiFi
// test BAM. Both publish to ${params.outdir}/fixtures so the test can compare the outputs against them.

// A reference with a second contig, chrUn_extra (a 1 Mb copy of chr19), and the BAM with the reads that fall wholly
// inside that window moved onto it. chrUn_extra plays the part of chrM or an unplaced contig: it has reads but no
// germline calls.
process PHASING_FIXTURE_BAM {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/fixtures", mode: 'copy'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)   // BAM aligned to fasta, which holds chr19 only

    output:
    tuple val(meta), path("fixture_input.bam"), path("fixture_input.bam.bai"), emit: bam
    tuple val(meta), path("reference.fa")                                     , emit: fasta
    tuple val(meta), path("reference.fa.fai")                                 , emit: fai

    script:
    def start = 54000001
    def end   = 55000000
    """
    zcat -f ${fasta} > chr19.fa
    samtools faidx chr19.fa
    samtools faidx chr19.fa chr19:${start}-${end} | sed '1s/.*/>chrUn_extra/' > extra.fa
    cat chr19.fa extra.fa > reference.fa
    rm chr19.fa chr19.fa.fai extra.fa
    samtools faidx reference.fa

    # Move every record whose alignment lies wholly inside the window onto chrUn_extra (SA tags dropped, as they
    # point at chr19 coordinates); everything else stays on chr19
    samtools view -h ${bam} \\
        | awk -F'\\t' -v OFS='\\t' -v s=${start} -v e=${end} '
            /^@SQ/ { print; if (\$2 == "SN:chr19") print "@SQ", "SN:chrUn_extra", "LN:" (e - s + 1); next }
            /^@/ { print; next }
            {
                span = 0; cig = \$6
                while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                    op = substr(cig, RLENGTH, 1)
                    if (op ~ /[MDN=X]/) span += substr(cig, 1, RLENGTH - 1)
                    cig = substr(cig, RLENGTH + 1)
                }
                if (\$3 == "chr19" && \$4 >= s && \$4 + span - 1 <= e) {
                    \$3 = "chrUn_extra"; \$4 = \$4 - s + 1
                    line = \$1; for (i = 2; i <= NF; i++) if (\$i !~ /^SA:Z:/) line = line OFS \$i
                    print line; next
                }
                print
            }' \\
        | samtools sort -o fixture_input.bam -
    samtools index fixture_input.bam

    if [ "\$(samtools view -c fixture_input.bam chrUn_extra)" -lt 20 ]; then
        echo "fixture_input.bam: fewer than 20 reads on chrUn_extra" >&2
        exit 1
    fi
    """
}

// Germline calls on chr19 (bcftools, so the test needs no Clair model) with a header that lists chr19 only, as
// ClairS-TO writes it, and one heterozygous SNV of them as the somatic VCF.
process PHASING_FIXTURE_VCF {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/fixtures", mode: 'copy'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.24--h118bc1c_2' :
        'biocontainers/bcftools:1.24--h118bc1c_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("germline.vcf.gz"), path("germline.vcf.gz.tbi"), emit: germline
    tuple val(meta), path("somatic.vcf.gz") , path("somatic.vcf.gz.tbi") , emit: somatic

    script:
    """
    zcat -f ${fasta} > ref.fa
    bcftools mpileup --threads ${task.cpus} -f ref.fa -r chr19 -a FORMAT/AD,FORMAT/DP -Ou ${bam} \\
        | bcftools call -mv -Ov -o called.vcf
    if grep '^##contig' called.vcf | grep -qv 'ID=chr19,'; then
        echo "called.vcf: header lists contigs other than chr19" >&2
        exit 1
    fi
    bcftools view -Oz -o germline.vcf.gz --write-index=tbi called.vcf
    bcftools view -i 'TYPE="snp" && GT="het"' called.vcf \\
        | awk '/^#/ { print; next } !done { print; done = 1 }' \\
        | bcftools view -Oz -o somatic.vcf.gz --write-index=tbi
    rm ref.fa ref.fa.fai called.vcf
    """
}
