#!/usr/bin/env python3
"""Generate the synthetic haplotagged BAM and marker BED for the OWL tests.

Three (AC)n loci on a 6 kb chr1, 6 reads per haplotype per locus, all in phase block 1000:
  locus 1  hap 1 stable (20 units), hap 2 unstable (16-26 units via CIGAR I/D, CV well above 5%)
  locus 2  both haplotypes stable
  locus 3  both haplotypes stable
Expected owl score: 1 of 6 haplotypes high, %high 16.67.
Writes owl.sam, then: samtools sort -o owl.bam owl.sam && samtools index owl.bam && rm owl.sam
"""
import os
import random

F = os.path.dirname(os.path.abspath(__file__))
random.seed(1)
CHROM_LEN = 6000
UNITS = 20
LOCI = [1000, 3000, 5000]
FLANK = 150
HAP2_UNSTABLE = [16, 18, 20, 22, 24, 26]

ref = [random.choice("ACGT") for _ in range(CHROM_LEN)]
for start in LOCI:
    ref[start:start + 2 * UNITS] = list("AC" * UNITS)
ref = "".join(ref)


def read(name, start, units, hap):
    """A read spanning the locus at `start` whose repeat has `units` AC units."""
    rstart = start - FLANK
    left = ref[rstart:start]
    right = ref[start + 2 * UNITS:start + 2 * UNITS + FLANK]
    seq = left + "AC" * units + right
    diff = 2 * (units - UNITS)
    if diff > 0:
        cigar = f"{FLANK + 2 * UNITS}M{diff}I{FLANK}M"
    elif diff < 0:
        cigar = f"{FLANK + 2 * units}M{-diff}D{FLANK}M"
    else:
        cigar = f"{len(seq)}M"
    return "\t".join([name, "0", "chr1", str(rstart + 1), "60", cigar, "*", "0", "0", seq, "I" * len(seq),
                      f"HP:i:{hap}", "PS:i:1000", "RG:Z:test"])


lines = ["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:chr1\tLN:{CHROM_LEN}", "@RG\tID:test\tSM:test"]
for i, start in enumerate(LOCI, 1):
    for j in range(6):
        lines.append(read(f"l{i}_h1_{j}", start, UNITS, 1))
        lines.append(read(f"l{i}_h2_{j}", start, HAP2_UNSTABLE[j] if i == 1 else UNITS, 2))

with open(os.path.join(F, "owl.sam"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
with open(os.path.join(F, "owl_markers.bed"), "w") as fh:
    fh.write("".join(f"chr1\t{s}\t{s + 2 * UNITS}\t(AC)n\n" for s in LOCI))
print("wrote owl.sam", len(lines) - 3, "reads; owl_markers.bed", len(LOCI), "loci")
