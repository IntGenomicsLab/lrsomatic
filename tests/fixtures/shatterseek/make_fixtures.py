#!/usr/bin/env python3
"""Generate the synthetic ASCAT / Wakhan / Severus / SAVANA fixtures for the SHATTERSEEK real-data tests.

One invented chromothripsis event in every caller's own file layout (GRCh38 coordinates, sample `test`):
  chr19:20.1-29.4 Mb  12 interleaved SVs, three each of DEL (+-), DUP (-+), h2hINV (++) and t2tINV (--)
                      (SV k joins breakpoints k and k+12, so every SV overlaps every other without nesting)
  chr19 <-> chr7      5 translocations from inside the cluster to chr7:50.1-50.9 Mb
  chr19:20-30 Mb      copy number alternating 1/2 over ten 1 Mb segments (2 elsewhere)
  chr2:10.00-10.05 Mb one isolated deletion (a chromosome with an SV but no event)
Expected ShatterSeek result: chr19 high confidence (HC1, HC2), chr2 and chr7 none.
Also present so the conversion is exercised: an insertion and a non-PASS deletion (both dropped), a Wakhan segment
with a fractional haplotype state and fractional SAVANA copy numbers (rounded), and SAVANA translocations whose
_2 mate (on chr7) comes first in the file while BP_NOTATION is written in (_1, _2) order.
"""
import gzip, os

F = os.path.dirname(os.path.abspath(__file__))
LEN = {"chr2": 242193529, "chr7": 159345973, "chr19": 58617616}
ORDER = {"chr2": 2, "chr7": 7, "chr19": 19}


def w(name, lines, gz=False):
    path = os.path.join(F, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with (gzip.open(path, "wt") if gz else open(path, "w")) as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote", name, len(lines), "lines")


def alt(strands, chrom, pos):
    """VCF breakend ALT for a record whose own strand is strands[0] and its mate's strands[1]."""
    own, mate = strands
    bracket = "[" if mate == "-" else "]"
    return f"N{bracket}{chrom}:{pos}{bracket}" if own == "+" else f"{bracket}{chrom}:{pos}{bracket}N"


# ---------------------------------------------------------------- the event
bps = [20_100_000 + k * 400_000 + (k % 3) * 7_000 for k in range(24)]
types = ["DEL", "DUP", "h2hINV", "t2tINV"] * 3
STRANDS = {"DEL": "+-", "DUP": "-+", "h2hINV": "++", "t2tINV": "--"}
intra = [("chr19", bps[k], bps[k + 12], types[k]) for k in range(12)]
tra = [(p19, p7, s) for p19, p7, s in zip([20_500_000, 22_500_000, 24_500_000, 26_500_000, 28_500_000],
                                          [50_100_000, 50_300_000, 50_500_000, 50_700_000, 50_900_000],
                                          ["+-", "-+", "++", "--", "+-"])]

# copy number (start, end, total): ten alternating 1 Mb segments inside chr19:20-30 Mb
cn = [("chr19", 1, 20_000_000, 2)]
cn += [("chr19", 20_000_001 + i * 1_000_000, 20_000_000 + (i + 1) * 1_000_000, 1 if i % 2 == 0 else 2) for i in range(10)]
cn += [("chr19", 30_000_001, LEN["chr19"], 2), ("chr2", 1, LEN["chr2"], 2), ("chr7", 1, LEN["chr7"], 2)]
cn.sort(key=lambda r: (ORDER[r[0]], r[1]))

# ---------------------------------------------------------------- Severus (somatic SV VCF)
sev_hdr = ["##fileformat=VCFv4.2", "##source=Severus_v1.6",
           '##FILTER=<ID=PASS,Description="All filters passed">',
           '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
           '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
           '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">',
           '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness">',
           '##INFO=<ID=MATE_ID,Number=1,Type=String,Description="MATE ID for breakends">',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">'] + \
          [f"##contig=<ID={c},length={n}>" for c, n in LEN.items()] + \
          ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_tumor"]
rec = []
def sev(chrom, pos, vid, alt_, info, flt="PASS"):
    rec.append((chrom, pos, f"{chrom}\t{pos}\t{vid}\tN\t{alt_}\t60.0\t{flt}\tPRECISE;{info}\tGT:DV\t0/1:12"))
for i, (c, p1, p2, t) in enumerate(intra, 1):
    s = STRANDS[t]
    if t in ("DEL", "DUP"):    # Severus writes deletions and duplications as symbolic alleles
        sev(c, p1, f"severus_{t}{i}", f"<{t}>", f"SVTYPE={t};SVLEN={p2 - p1};END={p2};STRANDS={s}")
    else:                      # and inversion-like junctions as mirrored breakend pairs
        sev(c, p1, f"severus_BND{i}_1", alt(s, c, p2), f"SVTYPE=BND;SVLEN={p2 - p1};MATE_ID=severus_BND{i}_2;STRANDS={s}")
        sev(c, p2, f"severus_BND{i}_2", alt(s[::-1], c, p1), f"SVTYPE=BND;SVLEN={p2 - p1};MATE_ID=severus_BND{i}_1;STRANDS={s[::-1]}")
for j, (p19, p7, s) in enumerate(tra, 1):
    sev("chr19", p19, f"severus_BND{100 + j}_1", alt(s, "chr7", p7), f"SVTYPE=BND;SVLEN=0;MATE_ID=severus_BND{100 + j}_2;STRANDS={s}")
    sev("chr7", p7, f"severus_BND{100 + j}_2", alt(s[::-1], "chr19", p19), f"SVTYPE=BND;SVLEN=0;MATE_ID=severus_BND{100 + j}_1;STRANDS={s[::-1]}")
sev("chr2", 10_000_000, "severus_DEL200", "<DEL>", "SVTYPE=DEL;SVLEN=50000;END=10050000;STRANDS=+-")
sev("chr2", 20_000_000, "severus_DEL201", "<DEL>", "SVTYPE=DEL;SVLEN=80000;END=20080000;STRANDS=+-", flt="FAIL_LOWSUPPORT")
sev("chr19", 25_000_100, "severus_INS1", "<INS>", "SVTYPE=INS;SVLEN=300")
rec.sort(key=lambda r: (ORDER[r[0]], r[1]))
w("severus/severus_somatic.vcf.gz", sev_hdr + [r[2] for r in rec], gz=True)

# ---------------------------------------------------------------- ASCAT (sample test.tumour; chromosomes without chr)
w("ascat/test.segments.txt", ["sample\tchr\tstartpos\tendpos\tnMajor\tnMinor"] +
  [f"test.tumour\t{c[3:]}\t{s}\t{e}\t1\t{t - 1}" for c, s, e, t in cn])

# ---------------------------------------------------------------- Wakhan 0.5.0 integer profile (top-ranked solution)
bed_hdr = ["#chr: chromosome number", "#start: start address for CN segment", "#end: end address for CN segment",
           "#hp1_coverage / hp2_coverage: median coverage for this segment, per haplotype",
           "#hp1_copynumber_state / hp2_copynumber_state: detected copy number state (integer/fraction), per haplotype",
           "#hp1_confidence / hp2_confidence: confidence score, per haplotype",
           "#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file (union across haplotypes)",
           "#chr\tstart\tend\thp1_coverage\thp1_copynumber_state\thp1_confidence\thp2_coverage\thp2_copynumber_state\thp2_confidence\tsvs_breakpoints_ids"]
def hp(c, t):
    return (1.0, 1.2) if c == "chr7" else (1.0, float(t - 1))   # chr7: a subclonal 1.2 state (total 2.2 -> 2)
w("wakhan/solution_rank_1/integer_profile.bed", bed_hdr +
  [f"{c}\t{s - 1}\t{e}\t{15 * hp(c, t)[0]:.2f}\t{hp(c, t)[0]}\t0.9\t{15 * hp(c, t)[1]:.2f}\t{hp(c, t)[1]}\t0.9\t" for c, s, e, t in cn])

# ---------------------------------------------------------------- SAVANA (prefix test_tumor; BND pairs, notation in (_1, _2) order)
sav_hdr = ["##fileformat=VCFv4.2", "##source=SAVANAv1.3.8",
           '##FILTER=<ID=PASS,Description="All filters passed">',
           '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
           '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">',
           '##INFO=<ID=TUMOUR_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumour reads">',
           '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
           '##INFO=<ID=BP_NOTATION,Number=1,Type=String,Description="+- notation format of variant (same for paired breakpoints)">',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'] + \
          [f"##contig=<ID={c},length={n}>" for c, n in LEN.items()] + \
          ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_tumor"]
rec = []
def sav_pair(n, c1, p1, c2, p2, note, svlen):
    info = f"SVTYPE=BND;TUMOUR_READ_SUPPORT=12;SVLEN={svlen};BP_NOTATION={note}"
    rec.append((c1, p1, f"{c1}\t{p1}\tID_{n}_1\tN\t{alt(note, c2, p2)}\t.\tPASS\tMATEID=ID_{n}_2;{info}\tGT\t0/1"))
    rec.append((c2, p2, f"{c2}\t{p2}\tID_{n}_2\tN\t{alt(note[::-1], c1, p1)}\t.\tPASS\tMATEID=ID_{n}_1;{info}\tGT\t0/1"))
for i, (c, p1, p2, t) in enumerate(intra, 1):
    sav_pair(i, c, p1, c, p2, STRANDS[t], p2 - p1)
for j, (p19, p7, s) in enumerate(tra, 1):   # _1 on chr19, so the chr7 _2 record is written first
    sav_pair(100 + j, "chr19", p19, "chr7", p7, s, 0)
sav_pair(200, "chr2", 10_000_000, "chr2", 10_050_000, "+-", 50000)
rec.append(("chr19", 25_000_100, "chr19\t25000100\tID_300_1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;TUMOUR_READ_SUPPORT=12;SVLEN=300;BP_NOTATION=<INS>\tGT\t0/1"))
rec.sort(key=lambda r: (ORDER[r[0]], r[1]))
w("savana/test_tumor.classified.somatic.vcf", sav_hdr + [r[2] for r in rec])
frac = {1: 1.08, 2: 1.94}   # SAVANA's absolute copy number is not an integer
w("savana/test_tumor_segmented_absolute_copy_number.tsv",
  ["chromosome\tstart\tend\tsegment_id\tbin_count\tsum_of_bin_lengths\tweight\tcopyNumber\tminorAlleleCopyNumber\tmeanBAF\tno_hetSNPs"] +
  [f"{c}\t{s}\t{e}\t{c}_seg{i + 1}\t10\t{e - s + 1}\t10.0\t{frac[t]}\t{0.0 if t == 1 else 0.97}\t0.5\t100" for i, (c, s, e, t) in enumerate(cn)])
print("done")
