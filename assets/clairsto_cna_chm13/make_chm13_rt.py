#!/usr/bin/env python
"""Build a CHM13 replication-timing (RT) file for ClairS-TO Verdict from the hg38 one.

RT is a smooth signal, so every hg38 RT row is lifted to CHM13 with the UCSC hg38->chm13v2
chain and each CHM13 ASCAT locus takes the values of the nearest lifted RT point on the same
chromosome. Output has the exact layout of RT_G1000_hg38.txt:
    <tab>Chr<tab>Position<tab><15 cell lines>
    1_809641<tab>1<tab>809641<tab>...

usage: make_chm13_rt.py RT_G1000_hg38.txt hg38-chm13v2.over.chain.gz loci_dir loci_glob out.txt
"""
import glob
import gzip
import os
import sys

import numpy as np

rt_fn, chain_fn, loci_dir, loci_glob, out_fn = sys.argv[1:6]
MAX_DIST = 2_000_000  # beyond this, fill with the chromosome median instead of the neighbour

# ---- 1. parse chain into ungapped blocks per target (hg38) chromosome -----------------------
blocks = {}  # tName -> list of (tStart, tEnd, qName, qStart, qStrand, qSize)
opener = gzip.open if chain_fn.endswith(".gz") else open
with opener(chain_fn, "rt") as fh:
    cur = None
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith("chain"):
            f = line.split()
            cur = dict(tName=f[2], tStart=int(f[5]), qName=f[7], qSize=int(f[8]), qStrand=f[9], qStart=int(f[10]))
            t, q = cur["tStart"], cur["qStart"]
            continue
        f = line.split()
        size = int(f[0])
        blocks.setdefault(cur["tName"], []).append((t, t + size, cur["qName"], q, cur["qStrand"], cur["qSize"]))
        if len(f) == 3:
            t += size + int(f[1])
            q += size + int(f[2])

chain_idx = {}
for chrom, bl in blocks.items():
    bl.sort()
    ts = np.array([b[0] for b in bl], dtype=np.int64)
    te = np.array([b[1] for b in bl], dtype=np.int64)
    chain_idx[chrom] = (ts, te, bl)
print(f"chain: {sum(len(v) for v in blocks.values())} ungapped blocks on {len(blocks)} hg38 contigs", file=sys.stderr)


def lift(chrom, pos1):
    """1-based hg38 position -> (qName, 1-based CHM13 position) or None."""
    if chrom not in chain_idx:
        return None
    ts, te, bl = chain_idx[chrom]
    p = pos1 - 1
    i = np.searchsorted(ts, p, side="right") - 1
    if i < 0 or p >= te[i]:
        return None
    tS, _tE, qName, qS, strand, qSize = bl[i]
    off = p - tS
    if strand == "+":
        q = qS + off
    else:
        q = qSize - (qS + off) - 1
    return qName, q + 1


# ---- 2. lift the hg38 RT rows -----------------------------------------------------------------
lifted = {}  # CHM13 chrom -> (positions list, values list)
header = None
n_in = n_ok = 0
with open(rt_fn) as fh:
    for line in fh:
        if header is None:
            header = line.rstrip("\n")
            continue
        f = line.rstrip("\n").split("\t")
        n_in += 1
        chrom = f[1] if f[1].startswith("chr") else "chr" + f[1]
        res = lift(chrom, int(float(f[2])))  # R wrote some positions as 8e+06
        if res is None:
            continue
        n_ok += 1
        d = lifted.setdefault(res[0], ([], []))
        d[0].append(res[1])
        d[1].append(f[3:])
print(f"RT rows: {n_in}, lifted: {n_ok}", file=sys.stderr)

for chrom in list(lifted):
    pos = np.array(lifted[chrom][0], dtype=np.int64)
    val = np.array(lifted[chrom][1], dtype=float)
    order = np.argsort(pos)
    lifted[chrom] = (pos[order], val[order])

# ---- 3. assign nearest lifted RT to every CHM13 locus ----------------------------------------
n_out = n_far = 0
with open(out_fn, "w") as out:
    out.write(header + "\n")
    for fn in sorted(glob.glob(os.path.join(loci_dir, loci_glob))):
        loci = np.loadtxt(fn, dtype=str, usecols=(0, 1), ndmin=2)
        if loci.size == 0:
            continue
        chrom = loci[0, 0]
        chrom_num = chrom.replace("chr", "")
        if chrom_num == "Y" or chrom not in lifted:
            print(f"skip {fn}: no lifted RT for {chrom}", file=sys.stderr)
            continue
        lpos = loci[:, 1].astype(float).astype(np.int64)
        pos, val = lifted[chrom]
        med = np.median(val, axis=0)
        i = np.searchsorted(pos, lpos)
        i_lo = np.clip(i - 1, 0, len(pos) - 1)
        i_hi = np.clip(i, 0, len(pos) - 1)
        use_hi = np.abs(pos[i_hi] - lpos) < np.abs(pos[i_lo] - lpos)
        nearest = np.where(use_hi, i_hi, i_lo)
        dist = np.abs(pos[nearest] - lpos)
        vals = val[nearest]
        far = dist > MAX_DIST
        vals[far] = med
        n_far += int(far.sum())
        for p, v in zip(lpos, vals):
            out.write(f"{chrom_num}_{p}\t{chrom_num}\t{p}\t" + "\t".join(f"{x:.6f}" for x in v) + "\n")
            n_out += 1
        print(f"{chrom}: {len(lpos)} loci, {int(far.sum())} beyond {MAX_DIST} bp from a lifted RT point", file=sys.stderr)
print(f"written {n_out} CHM13 RT rows ({n_far} filled with chromosome medians) to {out_fn}", file=sys.stderr)
