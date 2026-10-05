# Synthetic SV / copy-number fixtures for the PADFOOT and RECONPLOT tests

Small, invented caller outputs that let the two modules run for real in CI (the pipeline test profiles produce no
copy-number caller output, so no caller pair forms there). Nothing here comes from a patient: the events are made up
and only the reference slice is real sequence.

| file                                                                                                                           | mimics                                                                                                       | read by                                              |
| ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------- |
| `chr19_1-2Mb.fasta` (+ `.fai`)                                                                                                 | GRCh38 chr19:1-2,000,000, record renamed `chr19` so coordinates stay genomic                                 | PADFOOT (`--ref`)                                    |
| `severus_somatic.vcf.gz`                                                                                                       | Severus 1.6 somatic VCF (`STRANDS`, `DETAILED_TYPE`, `CLUSTERID`, `MATE_ID`, `DV`/`VAF`)                     | both                                                 |
| `wakhan/solution_rank_1/integer_profile.vcf`                                                                                   | Wakhan 0.5.0 top-ranked integer profile VCF (altered segments only, `CN1`/`CN2`/`COV1`) | PADFOOT                                              |
| `wakhan/solution_rank_1/integer_profile.bed`, `wakhan/solutions_ranks.tsv`                                                     | Wakhan 0.5.0 integer profile BED (every segment, `hp1_`/`hp2_copynumber_state`) + ranks                      | RECONPLOT (BED) / PADFOOT (ranks as the ploidy file) |
| `test.segments.txt`, `test.purityploidy.txt`, `test.tumour_tumourBAF.txt`                                                      | ASCAT (chromosomes without `chr`)                                                                            | RECONPLOT                                            |
| `test_tumor.classified.somatic.vcf` / `.bedpe`                                                                                 | SAVANA 1.3.8 classified somatic calls (`BP_NOTATION`, `TUMOUR_AF`, `TUMOUR_ALT_HP`, `MATEID`)                | PADFOOT / RECONPLOT                                  |
| `test_tumor_segmented_absolute_copy_number.tsv`, `test_tumor_fitted_purity_ploidy.tsv`, `test_tumor_allele_counts_hetSNPs.bed` | SAVANA copy number, fit and het-SNP BAF                                                                      | both / RECONPLOT                                     |

Events (all chr19): a deletion over STK11 exons 4-10 with a matching one-copy loss in every CN caller, a tandem
duplication of chr19:1.56-1.70 Mb with a matching gain, a 300 bp templated insertion (copied from a repeat-free STK11
intron) and a 300 bp Alu insertion (copied from chr19:81,234-81,533, so RepeatMasker has something to classify), and
two breakend pairs (`++` and `+-`). Sample name `test`.

Regenerate with `python3 make_fixtures.py` (the reference slice itself comes from
`samtools faidx Homo_sapiens_assembly38.fasta chr19:1-2000000 | sed '1s/.*/>chr19/'`, then `samtools faidx`;
the Severus VCF is `bgzip`ped afterwards).

`custom_annotation/` holds the same chr19 slice's genes and repeats in the formats a user brings for a genome without bundled
Padfoot annotations (the `--padfoot_gff` / `--padfoot_rm` route, e.g. CHM13): `chr19_1-2Mb.gff3.gz`, a GENCODE-style GFF3 of
every gene in the slice re-expanded from Padfoot's bundled hg38 table (plus CDS rows, a shorter decoy transcript per
multi-exon gene and a lncRNA that Padfoot must ignore), and `chr19_1-2Mb_rm.fa.out.gz`, the slice's bundled repeats in
RepeatMasker `.out` layout. Annotating the fixture SVs with them gives the same genes and repeat classes as the bundled hg38
annotations.
