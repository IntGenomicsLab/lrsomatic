# Synthetic SV / copy-number fixtures for the PADFOOT and RECONPLOT tests

Small, invented caller outputs that let the two modules run for real in CI (the pipeline test profiles produce no
copy-number caller output, so no caller pair forms there). Nothing here comes from a patient: the events are made up
and only the reference slice is real sequence.

| file | mimics | read by |
| --- | --- | --- |
| `chr19_1-2Mb.fasta` (+ `.fai`) | GRCh38 chr19:1-2,000,000, record renamed `chr19` so coordinates stay genomic | PADFOOT (`--ref`) |
| `severus_somatic.vcf.gz` | Severus 1.6 somatic VCF (`STRANDS`, `DETAILED_TYPE`, `CLUSTERID`, `MATE_ID`, `DV`/`VAF`) | both |
| `test_2.0_0.9_0.9_wakhan_cna_integers.vcf` | Wakhan `solution_1/vcf_output/*_wakhan_cna_integers.vcf` (altered segments only, `CN1`/`CN2`/`COV1`) | PADFOOT |
| `test_2.0_0.9_0.9_copynumbers_segments_HP_{1,2}.bed`, `solutions_ranks.tsv` | Wakhan `bed_output/` + ranks | RECONPLOT |
| `test.segments.txt`, `test.purityploidy.txt`, `test.tumour_tumourBAF.txt` | ASCAT (chromosomes without `chr`) | RECONPLOT |
| `test_tumor.classified.somatic.vcf` / `.bedpe` | SAVANA 1.3.8 classified somatic calls (`BP_NOTATION`, `TUMOUR_AF`, `TUMOUR_ALT_HP`, `MATEID`) | PADFOOT / RECONPLOT |
| `test_tumor_segmented_absolute_copy_number.tsv`, `test_tumor_fitted_purity_ploidy.tsv`, `test_tumor_allele_counts_hetSNPs.bed` | SAVANA copy number, fit and het-SNP BAF | both / RECONPLOT |

Events (all chr19): a deletion over STK11 exons 4-10 with a matching one-copy loss in every CN caller, a tandem
duplication of chr19:1.56-1.70 Mb with a matching gain, a 300 bp templated insertion (copied from a repeat-free STK11
intron) and a 300 bp Alu insertion (copied from chr19:81,234-81,533, so RepeatMasker has something to classify), and
two breakend pairs (`++` and `+-`). Sample name `test`.

Regenerate with `python3 make_fixtures.py` (the reference slice itself comes from
`samtools faidx Homo_sapiens_assembly38.fasta chr19:1-2000000 | sed '1s/.*/>chr19/'`, then `samtools faidx`;
the Severus VCF is `bgzip`ped afterwards).
