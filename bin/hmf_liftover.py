#!/usr/bin/env python3
"""Lift the coordinate-bearing HMF (oncoanalyser) report resources from GRCh38 to CHM13.

The report's HMF section reads five resources; four carry GRCh38 coordinates:
  - Ensembl data cache (ensembl_gene_data.csv, ensembl_trans_exon_data.csv): gene, transcript,
    exon and coding bounds. Only canonical transcripts are kept (all the report reads), and one
    is kept only if EVERY exon lifts whole, onto one contig and chain, in rank order; otherwise
    the transcript is dropped (and the gene, with no transcript left). IDs, phases and ranks
    are build-independent and copied verbatim.
  - KnownHotspots VCF: POS lifted, REF re-checked against the target FASTA (mismatches dropped),
    ##contig rewritten from the FASTA index, ##assembly added; bgzipped.
  - ActionableCodingPanel BED: intervals lifted, names kept.
  - known_fusion_data CSV: only the Overrides regions IG_RANGE / ENHANCER_RANGE /
    THREE_PRIME_RANGE (orient;chrom;start;end); the orientation flips on an inverted chain.
The driver panel has no coordinates and is copied. A one-line REFERENCE file in the Ensembl
cache directory declares the target build, which is what the report checks, and
hmf_liftover_report.tsv counts what was lifted and dropped per resource.

The chain is read directly (UCSC chain format, target = GRCh38 = the chain's "t" side), so the
only dependency beyond the standard library is pysam, for the FASTA and bgzip.
"""

import argparse
import bisect
import collections
import csv
import gzip
import os
import shutil
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
REGION_KEYS = ("IG_RANGE", "ENHANCER_RANGE", "THREE_PRIME_RANGE")
# Allowed difference between an interval's lifted and original length (indels between builds)
LENGTH_SLACK_BP = 50
LENGTH_SLACK_FRAC = 0.1
# Long ranges (IG loci) are lifted by sampling, every this many bp
SAMPLE_STEP_BP = 5000


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def chr_name(c):
    return c if c.startswith("chr") else "chr" + c


class ChainMap:
    """Ungapped alignment blocks of a chain file, indexed by source contig."""

    def __init__(self, path):
        self.blocks = collections.defaultdict(list)  # tname -> [(tstart, tend, qstart, chain)]
        self.chains = []  # (score, qname, qsize, qstrand)
        with open_text(path) as fh:
            chain = None
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    score, tname, _tsize, tstrand, tstart = float(f[1]), f[2], f[3], f[4], int(f[5])
                    qname, qsize, qstrand, qstart = f[7], int(f[8]), f[9], int(f[10])
                    if tstrand != "+":
                        sys.exit(f"chain {f[12]}: source strand must be '+'")
                    chain = len(self.chains)
                    self.chains.append((score, qname, qsize, qstrand))
                    t, q = tstart, qstart
                    continue
                size = int(f[0])
                self.blocks[tname].append((t, t + size, q, chain))
                if len(f) == 3:
                    t += size + int(f[1])
                    q += size + int(f[2])
        for tname in self.blocks:
            self.blocks[tname].sort()
        self.starts = {t: [b[0] for b in bl] for t, bl in self.blocks.items()}

    def map_pos(self, chrom, pos):
        """0-based source position -> (contig, 0-based position, strand, chain) or None.
        Overlapping chains: the highest-scoring one wins."""
        chrom = chr_name(chrom)
        bl = self.blocks.get(chrom)
        if not bl:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos) - 1
        best = None
        # Blocks of different chains can overlap; a short look-back finds them
        for j in range(i, max(-1, i - 64), -1):
            t0, t1, q0, ch = bl[j]
            if t0 <= pos < t1 and (best is None or self.chains[ch][0] > self.chains[best[3]][0]):
                best = (t0, t1, q0, ch)
        if best is None:
            return None
        t0, _t1, q0, ch = best
        _score, qname, qsize, qstrand = self.chains[ch]
        q = q0 + (pos - t0)
        if qstrand == "-":
            q = qsize - 1 - q
        return qname, q, qstrand, ch

    def lift_interval(self, chrom, start, end):
        """Half-open [start, end) -> (contig, start, end, strand) when both ends map through one
        chain and the length is kept (within slack); else None."""
        a, b = self.map_pos(chrom, start), self.map_pos(chrom, end - 1)
        if a is None or b is None or a[3] != b[3]:
            return None
        lo, hi = min(a[1], b[1]), max(a[1], b[1]) + 1
        if abs((hi - lo) - (end - start)) > max(LENGTH_SLACK_BP, LENGTH_SLACK_FRAC * (end - start)):
            return None
        return a[0], lo, hi, a[2]

    def lift_range(self, chrom, start, end):
        """A long region (an IG locus): sample it, keep the samples on the majority contig and
        strand, span them. None when nothing maps."""
        hits = [self.map_pos(chrom, p) for p in list(range(start, end, SAMPLE_STEP_BP)) + [end - 1]]
        hits = [h for h in hits if h is not None]
        if not hits:
            return None
        (qname, strand), _n = collections.Counter((h[0], h[2]) for h in hits).most_common(1)[0]
        qs = [h[1] for h in hits if h[0] == qname and h[2] == strand]
        return qname, min(qs), max(qs) + 1, strand


class Report:
    def __init__(self):
        self.rows = []

    def add(self, resource, item, n):
        self.rows.append((resource, item, n))

    def write(self, path):
        with open(path, "w") as fh:
            fh.write("resource\titem\tcount\n")
            for r in self.rows:
                fh.write("\t".join(map(str, r)) + "\n")


def lift_ensembl(cm, src, dst, report):
    """Canonical transcripts lifted whole; genes re-spanned; REFERENCE written by the caller."""
    os.makedirs(dst, exist_ok=True)
    with open(os.path.join(src, "ensembl_trans_exon_data.csv"), newline="") as fh:
        rd = csv.DictReader(fh)
        ex_fields = rd.fieldnames
        by_tx = collections.defaultdict(list)
        n_rows = 0
        for r in rd:
            n_rows += 1
            if r["CanonicalTranscriptId"] == r["TransId"]:
                by_tx[r["TransId"]].append(r)
    with open(os.path.join(src, "ensembl_gene_data.csv"), newline="") as fh:
        rd = csv.DictReader(fh)
        gene_fields = rd.fieldnames
        genes = {r["GeneId"]: r for r in rd}

    kept_rows, gene_span, dropped = [], {}, collections.Counter()
    for _tid, rows in by_tx.items():
        gene = genes.get(rows[0]["GeneId"])
        if gene is None:
            dropped["no gene row"] += 1
            continue
        plain = not gene["Chromosome"].startswith("chr")
        lifted = []
        for r in rows:
            m = cm.lift_interval(gene["Chromosome"], int(r["ExonStart"]) - 1, int(r["ExonEnd"]))
            if m is None:
                break
            lifted.append((int(r["ExonRank"]), m))
        if len(lifted) != len(rows):
            dropped["exon did not lift"] += 1
            continue
        contigs = {m[0] for _rk, m in lifted}
        strands = {m[3] for _rk, m in lifted}
        if len(contigs) != 1 or len(strands) != 1:
            dropped["exons on several contigs or strands"] += 1
            continue
        flip = strands == {"-"}
        old_strand = int(rows[0]["Strand"])
        new_strand = -old_strand if flip else old_strand
        # Rank order along the new strand: ascending coordinates on +, descending on -
        order = [m[1] for _rk, m in sorted(lifted)]
        if order != sorted(order, reverse=(new_strand == -1)):
            dropped["exon order broken"] += 1
            continue
        qname = contigs.pop()
        coding = None
        if rows[0]["CodingStart"] not in ("NULL", ""):
            a = cm.map_pos(gene["Chromosome"], int(rows[0]["CodingStart"]) - 1)
            b = cm.map_pos(gene["Chromosome"], int(rows[0]["CodingEnd"]) - 1)
            if a is None or b is None or a[0] != qname or b[0] != qname:
                dropped["coding bounds did not lift"] += 1
                continue
            coding = (min(a[1], b[1]) + 1, max(a[1], b[1]) + 1)
        t_start = min(m[1] for _rk, m in lifted) + 1
        t_end = max(m[2] for _rk, m in lifted)
        by_rank = dict(lifted)
        for r in rows:
            m = by_rank[int(r["ExonRank"])]
            r = dict(r)
            r.update(Strand=str(new_strand), TransStart=str(t_start), TransEnd=str(t_end),
                     ExonStart=str(m[1] + 1), ExonEnd=str(m[2]))
            if coding:
                r.update(CodingStart=str(coding[0]), CodingEnd=str(coding[1]))
            kept_rows.append(r)
        out_chrom = qname[3:] if plain and qname.startswith("chr") else qname
        gene_span[gene["GeneId"]] = (out_chrom, new_strand, t_start, t_end)

    with open(os.path.join(dst, "ensembl_trans_exon_data.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ex_fields, lineterminator="\n")
        w.writeheader()
        w.writerows(kept_rows)
    n_genes_dropped = 0
    with open(os.path.join(dst, "ensembl_gene_data.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=gene_fields, lineterminator="\n")
        w.writeheader()
        for gid, g in genes.items():
            if gid not in gene_span:
                n_genes_dropped += 1
                continue
            chrom, strand, start, end = gene_span[gid]
            # The gene spans its canonical transcript: the only part of it the report uses
            g = dict(g, Chromosome=chrom, Strand=str(strand), GeneStart=str(start), GeneEnd=str(end))
            w.writerow(g)
    report.add("ensembl", "canonical transcripts in", len(by_tx))
    report.add("ensembl", "canonical transcripts lifted", len(gene_span))
    for why, n in sorted(dropped.items()):
        report.add("ensembl", f"transcripts dropped: {why}", n)
    report.add("ensembl", "genes dropped (no lifted canonical transcript)", n_genes_dropped)
    report.add("ensembl", "non-canonical exon rows not carried", n_rows - sum(len(v) for v in by_tx.values()))
    return {g["GeneName"] for gid, g in genes.items() if gid not in gene_span}


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def lift_hotspots(cm, src, dst, fasta, assembly, report):
    fa = pysam.FastaFile(fasta)
    header, records, dropped = [], [], collections.Counter()
    with open_text(src) as fh:
        for line in fh:
            if line.startswith("##"):
                if not (line.startswith("##contig=") or line.startswith("##assembly=") or line.startswith("##reference=")):
                    header.append(line.rstrip("\n"))
                continue
            if line.startswith("#"):
                chrom_line = line.rstrip("\n")
                continue
            f = line.rstrip("\n").split("\t")
            pos0, ref, alts = int(f[1]) - 1, f[3], f[4].split(",")
            m = cm.lift_interval(f[0], pos0, pos0 + len(ref))
            if m is None:
                dropped["did not lift"] += 1
                continue
            qname, start, _end, strand = m
            if strand == "-":
                if len(ref) != 1 or any(len(a) != 1 for a in alts):
                    dropped["indel on an inverted chain"] += 1
                    continue
                ref, alts = revcomp(ref), [revcomp(a) for a in alts]
            if qname not in fa.references:
                dropped["contig not in FASTA"] += 1
                continue
            if fa.fetch(qname, start, start + len(ref)).upper() != ref.upper():
                dropped["REF differs on the target"] += 1
                continue
            records.append((qname, start + 1, "\t".join([qname, str(start + 1), f[2], ref, ",".join(alts)] + f[5:])))
    order = {c: i for i, c in enumerate(fa.references)}
    records.sort(key=lambda r: (order.get(r[0], len(order)), r[1]))
    plain = dst[:-3] if dst.endswith(".gz") else dst
    with open(plain, "w") as out:
        out.write(header[0] + "\n")  # ##fileformat
        out.write(f"##assembly={assembly}\n")
        for c, ln in zip(fa.references, fa.lengths):
            out.write(f"##contig=<ID={c},length={ln}>\n")
        for h in header[1:]:
            out.write(h + "\n")
        out.write(chrom_line + "\n")
        for r in records:
            out.write(r[2] + "\n")
    if dst.endswith(".gz"):
        pysam.tabix_compress(plain, dst, force=True)
        os.remove(plain)
    report.add("hotspots", "records lifted", len(records))
    for why, n in sorted(dropped.items()):
        report.add("hotspots", f"records dropped: {why}", n)


def lift_bed(cm, src, dst, report):
    rows, dropped = [], 0
    with open_text(src) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            m = cm.lift_interval(f[0], int(f[1]), int(f[2]))
            if m is None:
                dropped += 1
                continue
            rows.append([m[0], m[1], m[2]] + f[3:])
    rows.sort(key=lambda r: (r[0], r[1]))
    opener = gzip.open if dst.endswith(".gz") else open
    with opener(dst, "wt") as out:
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    report.add("actionable_bed", "regions lifted", len(rows))
    report.add("actionable_bed", "regions dropped", dropped)


def lift_fusions(cm, src, dst, report):
    n_lifted, n_dropped = 0, 0
    with open(src, newline="") as fh, open(dst, "w", newline="") as out:
        rd = csv.DictReader(fh)
        w = csv.DictWriter(out, fieldnames=rd.fieldnames, lineterminator="\n")
        w.writeheader()
        for r in rd:
            ov = (r.get("Overrides") or "").split()
            new = []
            for kv in ov:
                key, _, val = kv.partition("=")
                if key not in REGION_KEYS:
                    new.append(kv)
                    continue
                orient, chrom, start, end = val.split(";")
                m = cm.lift_range(chrom, int(start) - 1, int(end))
                if m is None:
                    n_dropped += 1
                    sys.stderr.write(f"{r['Type']} {r.get('FiveGene')}-{r.get('ThreeGene')}: {key} did not lift, dropped\n")
                    continue
                o = int(orient)
                if m[3] == "-" and o != 0:
                    o = -o
                new.append(f"{key}={o};{m[0]};{m[1] + 1};{m[2]}")
                n_lifted += 1
            r["Overrides"] = " ".join(new)
            w.writerow(r)
    report.add("known_fusions", "override regions lifted", n_lifted)
    report.add("known_fusions", "override regions dropped", n_dropped)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chain", required=True, help="GRCh38 -> target chain file (.chain or .chain.gz)")
    ap.add_argument("--fasta", required=True, help="Target FASTA (indexed), for the hotspot REF check")
    ap.add_argument("--outdir", required=True, help="Output bundle root (common/, dna/sv/, dna/variants/)")
    ap.add_argument("--driver-panel")
    ap.add_argument("--ensembl-data-dir")
    ap.add_argument("--known-fusions")
    ap.add_argument("--hotspots")
    ap.add_argument("--actionable-bed")
    ap.add_argument("--reference", default="t2t", help="Build written to the Ensembl cache's REFERENCE file")
    ap.add_argument("--assembly", default="T2T-CHM13v2.0", help="##assembly= of the lifted hotspot VCF")
    a = ap.parse_args()

    cm = ChainMap(a.chain)
    report = Report()
    common = os.path.join(a.outdir, "common")
    sv = os.path.join(a.outdir, "dna", "sv")
    var = os.path.join(a.outdir, "dna", "variants")
    for d in (common, sv, var):
        os.makedirs(d, exist_ok=True)

    if a.driver_panel:
        shutil.copyfile(a.driver_panel, os.path.join(common, os.path.basename(a.driver_panel)))
    if a.ensembl_data_dir:
        ens = os.path.join(common, "ensembl_data")
        lost = lift_ensembl(cm, a.ensembl_data_dir, ens, report)
        with open(os.path.join(ens, "REFERENCE"), "w") as fh:
            fh.write(a.reference + "\n")
        if a.driver_panel:
            with open(a.driver_panel) as fh:
                panel = {line.split("\t")[0] for line in list(fh)[1:]}
            dropped_panel = sorted(panel & lost)
            report.add("driver_panel", "panel genes without a lifted transcript", len(dropped_panel))
            for g in dropped_panel:
                report.add("driver_panel", f"not lifted: {g}", 1)
    if a.known_fusions:
        lift_fusions(cm, a.known_fusions, os.path.join(sv, os.path.basename(a.known_fusions)), report)
    if a.hotspots:
        name = os.path.basename(a.hotspots)
        lift_hotspots(cm, a.hotspots, os.path.join(var, name if name.endswith(".gz") else name + ".gz"),
                      a.fasta, a.assembly, report)
    if a.actionable_bed:
        lift_bed(cm, a.actionable_bed, os.path.join(var, os.path.basename(a.actionable_bed)), report)
    report.write(os.path.join(a.outdir, "hmf_liftover_report.tsv"))


if __name__ == "__main__":
    main()
