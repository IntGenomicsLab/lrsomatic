# CHM13 replication timing for ClairS-TO Verdict

ClairS-TO's Verdict module needs four resource files per reference build: per-chromosome 1000 Genomes
loci and allele files, a GC-content table and a replication-timing (RT) table. For CHM13 the pipeline
already downloads the first three as its ASCAT reference files (`ascat_loci`, `ascat_alleles`,
`ascat_loci_gc` in `conf/igenomes.config`) and `CLAIRSTO_CNA_RESOURCES` assembles them into the layout
ClairS-TO expects. The RT table is the only piece that has to be generated, which is what
`make_chm13_rt.py` does.

Upstream ships GRCh38 resources only, and the only RT track available is the one inside the ClairS-TO
image. RT is a smooth, large-scale signal, so lifting it is reasonable: every GRCh38 RT row is mapped
to CHM13 with the UCSC `hg38-chm13v2` chain, and each CHM13 ASCAT locus takes the values of the
nearest lifted point on the same chromosome (falling back to the chromosome median beyond 2 Mb).

## Regenerating `RT_CHM13.txt`

```bash
# 1. the GRCh38 RT table, from any ClairS-TO image
apptainer exec <clairs-to.sif> \
    cat /opt/micromamba/envs/clairs-to/bin/clairs-to_cna_data/reference_files/RT_G1000_hg38.txt \
    > RT_G1000_hg38.txt

# 2. lift it onto the CHM13 ASCAT loci
python make_chm13_rt.py \
    RT_G1000_hg38.txt \
    hg38-chm13v2.over.chain.gz \
    <dir with the unzipped G1000_loci_CHM13 files> \
    'G1000_loci_CHM13_chr*.txt' \
    RT_CHM13.txt

# 3. publish as the zip referenced by `clairsto_cna_rt` in conf/igenomes.config
zip RT_CHM13.txt.zip RT_CHM13.txt
```

The output must keep bare contig names in the `Chr` column (`1`, not `chr1`) — Verdict's
`correct_logr.py` prepends the `chr` prefix itself.

The chain file is the UCSC `hg38-chm13v2.over.chain.gz` from the T2T consortium.
