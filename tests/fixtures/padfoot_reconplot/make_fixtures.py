#!/usr/bin/env python3
"""Generate the synthetic Severus / Wakhan / ASCAT / SAVANA fixtures for the PADFOOT and RECONPLOT real-data tests.

Everything lives on GRCh38 chr19:1-2,000,000 (the reference slice chr19_1-2Mb.fasta, built with
`samtools faidx Homo_sapiens_assembly38.fasta chr19:1-2000000 | sed '1s/.*/>chr19/'`). The events are
invented but use the exact file layouts, column names and INFO/FORMAT keys the two modules read:
  DEL   chr19:1,215,001-1,223,000  STK11 exons 4-10, with a matching one-copy loss (LOH) in every CN caller
  DUP   chr19:1,560,001-1,700,000  MBD3/TCF3 region, tandem duplication with a matching 3-copy gain
  INS   chr19:1,100,500            300 bp copied from a repeat-free STK11 intron (templated insertion)
  INS   chr19:1,050,100            300 bp Alu copied from chr19:81,234-81,533 (RepeatMasker has something to find)
  BND   chr19:400,000 <-> 900,000  ++ breakend pair (inversion-like)
  BND   chr19:1,050,000 <-> 1,900,000  +- breakend pair
Sample name: `test` (pipeline meta.id); Severus/SAVANA tumour column `test_tumor`; ASCAT `test.tumour`;
Wakhan solution `2.0_0.9_0.9`. Chromosome length of chr19 in headers is the real one (58,617,616).
"""
import os, sys
F = os.path.dirname(os.path.abspath(__file__))
CHR, CHR_LEN = "chr19", 58617616

def read_slice(start, end):   # 1-based inclusive, from the bundled fasta slice
    seq = []
    with open(os.path.join(F, "chr19_1-2Mb.fasta")) as fh:
        next(fh)
        for line in fh:
            seq.append(line.strip())
    s = "".join(seq)
    return s[start - 1:end].upper()

def ref_base(pos):
    return read_slice(pos, pos)

def w(name, lines):
    with open(os.path.join(F, name), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote", name, len(lines), "lines")

ins_unique = read_slice(1181701, 1182000)   # repeat-free window in STK11 intron 3
ins_alu    = read_slice(81234, 81533)       # SINE/Alu element
assert len(ins_unique) == 300 and len(ins_alu) == 300 and "N" not in ins_unique + ins_alu

# ---------------------------------------------------------------- Severus (somatic SV VCF, pipeline name severus_somatic.vcf.gz)
sev_hdr = f"""##fileformat=VCFv4.2
##source=Severus_v1.6
##fileDate=2026-10-01
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Reciprocal Inversion">
##ALT=<ID=BND,Description="Breakend">
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="SV with precise breakpoints coordinates and length">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="SV with imprecise breakpoints coordinates and length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness">
##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description="Detailed type of the SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Length of the unmapped sequence between breakpoint">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of supporting reads">
##INFO=<ID=PHASESETID,Number=1,Type=String,Description="Matching phaseset ID for phased SVs">
##INFO=<ID=HP,Number=1,Type=Integer,Description="Matching haplotype ID for phased SVs">
##INFO=<ID=CLUSTERID,Number=1,Type=String,Description="Cluster ID in breakpoint_graph">
##INFO=<ID=INSSEQ,Number=1,Type=String,Description="Insertion sequence between breakpoints">
##INFO=<ID=MATE_ID,Number=1,Type=String,Description="MATE ID for breakends">
##INFO=<ID=INSIDE_VNTR,Number=1,Type=String,Description="True if an indel is inside a VNTR">
##INFO=<ID=ALIGNED_POS,Number=1,Type=String,Description="Position in the reference">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##FORMAT=<ID=hVAF,Number=3,Type=Float,Description="Haplotype specific variant Allele frequency (H0,H1,H2)">
##contig=<ID={CHR},length={CHR_LEN}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_tumor""".split("\n")
FMT = "GT:VAF:hVAF:DR:DV"
def sev(pos, vid, alt, info, vaf, dr, dv):
    return f"{CHR}\t{pos}\t{vid}\tN\t{alt}\t60.0\tPASS\t{info}\t{FMT}\t0/1:{vaf:.2f}:{vaf:.2f},0.00,0.00:{dr}:{dv}"
sev_rec = [
    sev(400000,  "severus_BND1_1", f"N]{CHR}:900000]",  "PRECISE;SVTYPE=BND;SVLEN=500000;MATE_ID=severus_BND1_2;STRANDS=++;DETAILED_TYPE=complex_inv;MAPQ=60.0;CLUSTERID=severus_3", 0.35, 26, 14),
    sev(900000,  "severus_BND1_2", f"N]{CHR}:400000]",  "PRECISE;SVTYPE=BND;SVLEN=500000;MATE_ID=severus_BND1_1;STRANDS=++;DETAILED_TYPE=complex_inv;MAPQ=60.0;CLUSTERID=severus_3", 0.35, 26, 14),
    sev(1050000, "severus_BND2_1", f"N[{CHR}:1900000[", "PRECISE;SVTYPE=BND;SVLEN=850000;MATE_ID=severus_BND2_2;STRANDS=+-;MAPQ=60.0;CLUSTERID=severus_4", 0.30, 28, 12),
    sev(1050100, "severus_INS2",   ins_alu,             "PRECISE;SVTYPE=INS;SVLEN=300;MAPQ=60.0", 0.38, 31, 19),
    sev(1100500, "severus_INS1",   ins_unique,          "PRECISE;SVTYPE=INS;SVLEN=300;MAPQ=60.0", 0.40, 30, 20),
    sev(1215001, "severus_DEL1",   "<DEL>",             "PRECISE;SVTYPE=DEL;SVLEN=8000;END=1223000;STRANDS=+-;MAPQ=60.0;CLUSTERID=severus_1", 0.45, 22, 18),
    sev(1560001, "severus_DUP1",   "<DUP>",             "PRECISE;SVTYPE=DUP;SVLEN=140000;END=1700000;STRANDS=-+;DETAILED_TYPE=tandem_duplication;MAPQ=60.0;CLUSTERID=severus_2", 0.60, 16, 24),
    sev(1900000, "severus_BND2_2", f"]{CHR}:1050000]N", "PRECISE;SVTYPE=BND;SVLEN=850000;MATE_ID=severus_BND2_1;STRANDS=+-;MAPQ=60.0;CLUSTERID=severus_4", 0.30, 28, 12),
]
w("severus_somatic.vcf", sev_hdr + sev_rec)   # bgzipped afterwards

# ---------------------------------------------------------------- Wakhan (solution 2.0_0.9_0.9)
SOL = "2.0_0.9_0.9"
wak_hdr = f"""##fileformat=VCFv4.2
##source=Wakhan_v0.4.2
##fileDate=2026-10-01
##ALT=<ID=CNV,Description="Copy number variant region">
##ALT=<ID=DEL,Description="Deletion relative to the reference">
##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">
##INFO=<ID=REFLEN,Number=1,Type=Integer,Description="Number of REF positions included in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=HET,Number=0,Type=Flag,Description="Segment is heterogeneous">
##INFO=<ID=BPS,Number=0,Type=String,Description="Breakpoints covering segment">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=TCN,Number=1,Type=Float,Description="Estimated total copy numbers of segment">
##FORMAT=<ID=CN1,Number=1,Type=Float,Description="Estimated haplotype-1 segment copy number">
##FORMAT=<ID=CN2,Number=1,Type=Float,Description="Estimated haplotype-2 segment copy number">
##FORMAT=<ID=CNQ1,Number=1,Type=Float,Description="Estimated haplotype-1 segment confidence score">
##FORMAT=<ID=CNQ2,Number=1,Type=Float,Description="Estimated haplotype-2 segment confidence score">
##FORMAT=<ID=COV1,Number=1,Type=Float,Description="Estimated haplotype-1 segment coverage value">
##FORMAT=<ID=COV2,Number=1,Type=Float,Description="Estimated haplotype-2 segment coverage value">
##contig=<ID={CHR},length={CHR_LEN}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample""".split("\n")
WF = "GT:TCN:CN1:CN2:CNQ1:CNQ2:COV1:COV2"
wak_rec = [
    f"{CHR}\t1215001\twakhan:LOSS:{CHR}:1215001-1223000\tN\t<DEL>\t1000\tPASS\tSVTYPE=CNV;SVLEN=7999;END=1223000;BPS=severus_DEL1\t{WF}\t1/1:1.0:1.0:0.0:0.91:0.91:15.0:0.0",
    f"{CHR}\t1560001\twakhan:GAIN:{CHR}:1560001-1700000\tN\t<DUP>\t1000\tPASS\tSVTYPE=CNV;SVLEN=139999;END=1700000;BPS=severus_DUP1\t{WF}\t0/1:3.0:2.0:1.0:0.88:0.90:30.0:15.0",
]
w(f"test_{SOL}_wakhan_cna_integers.vcf", wak_hdr + wak_rec)
w("solutions_ranks.tsv", ["repository_name\tdna_purity\tcell_purity\tploidy\tconfidence\tsolution_rank", f"{SOL}\t0.9\t0.9\t2.0\t0.9\t1"])
bed_hdr = ["#chr: chromosome number", "#start: start address for CN segment", "#end: end address for CN segment",
           "#coverage: median coverage for this segment", "#copynumber_state: detected copy number state (integer/fraction)",
           "#confidence: confidence score", "#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file",
           "#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids"]
segs = [(0, 1215000, 1, 1, "[]"), (1215001, 1223000, 1, 0, "['severus_DEL1']"), (1223001, 1560000, 1, 1, "[]"),
        (1560001, 1700000, 2, 1, "['severus_DUP1']"), (1700001, CHR_LEN, 1, 1, "[]")]
for hp, idx in (("HP_1", 2), ("HP_2", 3)):
    rows = [f"{CHR}\t{s}\t{e}\t{15.0 * seg[idx]:.2f}\t{float(seg[idx])}\t0.9\t{seg[4]}" for seg in segs for s, e in [(seg[0], seg[1])]]
    w(f"test_{SOL}_copynumbers_segments_{hp}.bed", bed_hdr + rows)

# ---------------------------------------------------------------- ASCAT (sample test.tumour; chromosomes without the chr prefix)
asc = [(1, 1215000, 1, 1), (1215001, 1223000, 1, 0), (1223001, 1560000, 1, 1), (1560001, 1700000, 2, 1), (1700001, CHR_LEN, 1, 1)]
w("test.segments.txt", ["sample\tchr\tstartpos\tendpos\tnMajor\tnMinor"] + [f"test.tumour\t19\t{s}\t{e}\t{a}\t{b}" for s, e, a, b in asc])
w("test.purityploidy.txt", ["AberrantCellFraction\tPloidy", "0.9\t2.0"])
def baf_at(pos):
    if 1215001 <= pos <= 1223000: return 1.0 if pos % 100000 < 50000 else 0.0
    if 1560001 <= pos <= 1700000: return 0.333 if pos % 100000 < 50000 else 0.667
    return 0.5
snps = list(range(25000, 2000000, 25000))
w("test.tumour_tumourBAF.txt", ["\tChromosome\tPosition\ttest.tumour"] + [f"19_{p}\t19\t{p}\t{baf_at(p)}" for p in snps])

# ---------------------------------------------------------------- SAVANA (prefix test_tumor)
sav_hdr = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20261001
##source=SAVANAv1.3.8
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">
##INFO=<ID=NORMAL_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting normal reads">
##INFO=<ID=TUMOUR_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumour reads">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=TUMOUR_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of tumour variant-supporting reads to tumour read depth (DP) at breakpoint">
##INFO=<ID=NORMAL_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of normal variant-supporting reads to normal read depth (DP) at breakpoint">
##INFO=<ID=BP_NOTATION,Number=1,Type=String,Description="+- notation format of variant (same for paired breakpoints)">
##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of evidence for a breakpoint - CIGAR (INS, DEL, SOFTCLIP), SUPPLEMENTARY or mixture">
##INFO=<ID=TUMOUR_ALT_HP,Number=3,Type=Integer,Description="Counts of SV-supporting reads belonging to each haplotype in the tumour sample (1/2/NA)">
##INFO=<ID=NORMAL_ALT_HP,Number=3,Type=Integer,Description="Counts of reads belonging to each haplotype in the normal sample (1/2/NA)">
##INFO=<ID=CLASS,Number=1,Type=String,Description="Variant class prediction from model ont-somatic.pkl">
##contig=<ID={CHR},length={CHR_LEN}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_tumor""".split("\n")
def sav(pos, vid, alt, info):
    return f"{CHR}\t{pos}\t{vid}\t{ref_base(pos)}\t{alt}\t.\tPASS\t{info}\tGT\t0/1"
def bnd(pos, vid, mate, alt_tpl, svlen, note, supp, af, hp):
    b = ref_base(pos)
    alt = alt_tpl.replace("B", b)
    return sav(pos, vid, alt, f"SVTYPE=BND;MATEID={mate};TUMOUR_READ_SUPPORT={supp};NORMAL_READ_SUPPORT=0;SVLEN={svlen};BP_NOTATION={note};SOURCE=SUPPLEMENTARY;TUMOUR_AF={af};NORMAL_AF=0,0;TUMOUR_ALT_HP={hp};NORMAL_ALT_HP=0,0,0;CLASS=PREDICTED_SOMATIC")
sav_rec = [
    bnd(400000,  "ID_4_1", "ID_4_2", f"B]{CHR}:900000]",  500000, "++", 14, "0.35,0.34", "8,2,4"),
    bnd(900000,  "ID_4_2", "ID_4_1", f"B]{CHR}:400000]",  500000, "++", 14, "0.34,0.35", "8,2,4"),
    bnd(1050000, "ID_5_1", "ID_5_2", f"[{CHR}:1900000[B", 850000, "--", 12, "0.30,0.29", "0,9,3"),
    sav(1100500, "ID_3_1", "<INS>", "SVTYPE=INS;TUMOUR_READ_SUPPORT=20;NORMAL_READ_SUPPORT=0;SVLEN=300;BP_NOTATION=<INS>;SOURCE=CIGAR;TUMOUR_AF=0.40,0.40;NORMAL_AF=0,0;TUMOUR_ALT_HP=0,0,20;NORMAL_ALT_HP=0,0,0;CLASS=PREDICTED_SOMATIC"),
    bnd(1215001, "ID_1_1", "ID_1_2", f"B[{CHR}:1223000[", 7999,   "+-", 18, "0.45,0.44", "12,2,4"),
    bnd(1223000, "ID_1_2", "ID_1_1", f"]{CHR}:1215001]B", 7999,   "+-", 18, "0.44,0.45", "12,2,4"),
    bnd(1560001, "ID_2_1", "ID_2_2", f"]{CHR}:1700000]B", 139999, "-+", 24, "0.60,0.58", "16,4,4"),
    bnd(1700000, "ID_2_2", "ID_2_1", f"B[{CHR}:1560001[", 139999, "-+", 24, "0.58,0.60", "16,4,4"),
    bnd(1900000, "ID_5_2", "ID_5_1", f"[{CHR}:1050000[B", 850000, "--", 12, "0.29,0.30", "0,9,3"),
]
w("test_tumor.classified.somatic.vcf", sav_hdr + sav_rec)
w("test_tumor.classified.somatic.bedpe", [
    f"{CHR}\t400000\t400000\t{CHR}\t900000\t900000\tID_4|500000bp|TUMOUR_14|++",
    f"{CHR}\t1050000\t1050000\t{CHR}\t1900000\t1900000\tID_5|850000bp|TUMOUR_12|--",
    f"{CHR}\t1100500\t1100500\t{CHR}\t1100501\t1100501\tID_3|300bp|TUMOUR_20|<INS>",
    f"{CHR}\t1215001\t1215001\t{CHR}\t1223000\t1223000\tID_1|7999bp|TUMOUR_18|+-",
    f"{CHR}\t1560001\t1560001\t{CHR}\t1700000\t1700000\tID_2|139999bp|TUMOUR_24|-+",
])
cn_rows = [(1, 1215000, 121, 2.0, 1.0, 0.5, 400), (1215001, 1223000, 1, 1.0, 0.0, 1.0, 5), (1223001, 1560000, 34, 2.0, 1.0, 0.5, 120),
           (1560001, 1700000, 14, 3.0, 1.0, 0.667, 50), (1700001, CHR_LEN, 5692, 2.0, 1.0, 0.5, 18000)]
w("test_tumor_segmented_absolute_copy_number.tsv",
  ["chromosome\tstart\tend\tsegment_id\tbin_count\tsum_of_bin_lengths\tweight\tcopyNumber\tminorAlleleCopyNumber\tmeanBAF\tno_hetSNPs"] +
  [f"{CHR}\t{s}\t{e}\t{CHR}_seg{i+1}\t{n}\t{e-s+1}\t{float(n)}\t{cn}\t{mn}\t{baf}\t{k}" for i, (s, e, n, cn, mn, baf, k) in enumerate(cn_rows)])
w("test_tumor_fitted_purity_ploidy.tsv", ["purity\tploidy\tdistance\trank", "0.9\t2.0\t0.1\t1"])
w("test_tumor_allele_counts_hetSNPs.bed", [f"{CHR}\t{p}\t{p}\tA\tG\t0\t10\t0\t10\t0\t{1 - baf_at(p):.3f}\t{baf_at(p):.3f}\t{CHR}_0" for p in snps])
print("done")
