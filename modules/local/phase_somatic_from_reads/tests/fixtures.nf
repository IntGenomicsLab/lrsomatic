// Test helper, not part of the pipeline: a 300 bp reference, haplotagged reads as BAM and CRAM, and a somatic VCF
// with one record per case of PHASE_SOMATIC_FROM_READS (see the table in main.nf.test).
process PSFR_FIXTURE {
    tag "fixture"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0':
        'biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("somatic.vcf.gz"), path("somatic.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("tumour.bam") , path("tumour.bam.bai")  , emit: bam
    tuple val(meta), path("tumour.cram"), path("tumour.cram.crai"), emit: cram
    tuple val(meta), path("ref.fa")                               , emit: fasta
    tuple val(meta), path("ref.fa.fai")                           , emit: fai

    script:
    """
    #!/usr/bin/env python3
    import random
    import pysam

    rng = random.Random(7)
    ref = ''.join(rng.choice('ACGT') for _ in range(300))
    with open('ref.fa', 'w') as fh:
        fh.write('>chrT\\n' + ref + '\\n')
    pysam.faidx('ref.fa')
    sub = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}

    def snv(seq, pos):
        return seq[:pos - 1] + sub[ref[pos - 1]] + seq[pos:]

    # Reads span the whole contig. HP1 reads r1-r4, HP2 reads r5-r8, untagged reads r9-r12; all tagged reads PS 1001.
    #   HP1: SNV 50, 230 and 250; deletion of 2 bp after 150; SNV 110 on r1-r3, SNV 90 on r1-r2
    #   HP2: SNV 70 and 250; insertion TT after 200 and GG after 230; SNV 110 on r5-r6
    #   untagged: SNV 130 and 250
    reads = []
    for i in range(1, 13):
        seq = ref
        hp = 1 if i <= 4 else 2 if i <= 8 else None
        if hp == 1:
            for p in (50, 230, 250):
                seq = snv(seq, p)
            if i <= 3:
                seq = snv(seq, 110)
            if i <= 2:
                seq = snv(seq, 90)
            seq = seq[:150] + seq[152:]
            cigar = '150M2D148M'
        elif hp == 2:
            for p in (70, 250):
                seq = snv(seq, p)
            if i <= 6:
                seq = snv(seq, 110)
            seq = seq[:200] + 'TT' + seq[200:230] + 'GG' + seq[230:]
            cigar = '200M2I30M2I70M'
        else:
            for p in (130, 250):
                seq = snv(seq, p)
            cigar = '300M'
        reads.append((f'r{i}', seq, cigar, hp))

    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chrT', 'LN': 300}],
              'RG': [{'ID': 'rg1', 'SM': 'tumour'}]}
    for name, mode in (('tumour.bam', 'wb'), ('tumour.cram', 'wc')):
        with pysam.AlignmentFile(name, mode, header=header, reference_filename='ref.fa') as out:
            for qname, seq, cigar, hp in reads:
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.reference_id, a.reference_start = qname, 0, 0
                a.mapping_quality, a.flag, a.cigarstring = 60, 0, cigar
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
                tags = [('RG', 'rg1')]
                if hp:
                    tags += [('HP', hp), ('PS', 1001)]
                a.set_tags(tags)
                out.write(a)
        pysam.index(name)

    a = lambda p: sub[ref[p - 1]]
    records = [
        (50, ref[49], a(50), '0/1'),
        (70, ref[69], a(70), '0/1'),
        (90, ref[89], a(90), '0/1'),
        (110, ref[109], a(110), '0/1'),
        (130, ref[129], a(130), '0/1'),
        (150, ref[149:152], ref[149], '0/1'),
        (200, ref[199], ref[199] + 'TT', '0/1'),
        (230, ref[229], a(230), '0/1'),
        (230, ref[229], ref[229] + 'GG', '0/1'),
        (250, ref[249], a(250), '1/1'),
        (270, ref[269], a(270), '0/0'),
    ]
    with open('somatic.vcf', 'w') as fh:
        fh.write('##fileformat=VCFv4.2\\n##contig=<ID=chrT,length=300>\\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n')
        fh.write('#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\ttumour\\n')
        for pos, r, alt, gt in records:
            fh.write(f'chrT\\t{pos}\\t.\\t{r}\\t{alt}\\t.\\tPASS\\t.\\tGT\\t{gt}\\n')
    pysam.tabix_index('somatic.vcf', preset='vcf', force=True)
    """
}
