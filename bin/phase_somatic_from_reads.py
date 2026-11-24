#!/usr/bin/env python3
"""Phase somatic small variants from the HP/PS tags of their alt-supporting reads in a haplotagged BAM/CRAM.

The haplotags come from the germline phase, so the somatic GT orientation and PS written here agree with the
phased germline VCF and with the reads by construction. Only alt reads vote: a somatic allele sits on one
haplotype, and ref reads come from both haplotypes (and from normal cells), so they carry no information about it.

Per biallelic SNV, deletion or insertion whose genotype carries the ALT allele:
  FORMAT/AHP = alt reads on HP1, on HP2, and without a usable HP/PS tag.
  A heterozygous record (0/1) is phased when its alt reads with a (PS, HP) tag number at least --min-reads and one
  (PS, HP) holds at least --min-fraction of them: GT 1|0 for HP1, 0|1 for HP2, PS = that read PS.
Every other record keeps its genotype, unphased and without PS. All records are written, in input order.
"""
import argparse
import multiprocessing
import sys

import pysam


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Somatic VCF (bgzipped or plain)")
    p.add_argument("--alignments", required=True, help="Haplotagged BAM/CRAM (HP and PS read tags), indexed")
    p.add_argument("--reference", required=True, help="Reference FASTA (needed to decode CRAM)")
    p.add_argument("--output", required=True, help="Output VCF, bgzipped (.vcf.gz); a .tbi index is written next to it")
    p.add_argument("--sample", help="Sample column to phase (default: the only sample in the VCF)")
    p.add_argument("--min-reads", type=int, default=3, help="Minimum HP/PS-tagged alt reads to phase a record")
    p.add_argument("--min-fraction", type=float, default=0.8, help="Minimum fraction of tagged alt reads on one (PS, HP)")
    p.add_argument("--min-mapq", type=int, default=10, help="Minimum mapping quality of a counted read")
    p.add_argument("--min-baseq", type=int, default=10, help="Minimum base quality of an SNV alt base")
    p.add_argument("--threads", type=int, default=1, help="Worker processes (one contig at a time each)")
    return p.parse_args(argv)


def variant_kind(ref, alt):
    """'snv', 'del', 'ins', or None for alleles this script does not match (MNVs, complex, symbolic)."""
    if alt is None or alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
        return None
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
        return "del"
    if len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
        return "ins"
    return None


def has_alt(gt):
    return gt is not None and any(a is not None and a > 0 for a in gt)


def is_het(gt):
    return gt is not None and len(gt) == 2 and set(gt) == {0, 1}


def read_hap(aln):
    """(PS, HP) of a haplotagged read, or None if it has no usable tag."""
    if not (aln.has_tag("HP") and aln.has_tag("PS")):
        return None
    hp = aln.get_tag("HP")
    if hp not in (1, 2):
        return None
    return aln.get_tag("PS"), hp


def alt_votes(bam, chrom, pos, ref, alt, kind, args):
    """{(PS, HP): alt reads} and the number of alt reads without a usable tag, at one record."""
    votes = {}
    untagged = 0
    start = pos - 1  # the VCF POS base (the anchor base of an indel), 0-based
    for col in bam.pileup(chrom, start, start + 1, truncate=True, stepper="nofilter", flag_filter=0,
                          min_base_quality=0, ignore_overlaps=False, ignore_orphans=False, max_depth=100000,
                          multiple_iterators=False):
        for pr in col.pileups:
            aln = pr.alignment
            if (aln.is_secondary or aln.is_supplementary or aln.is_duplicate or aln.is_qcfail
                    or aln.is_unmapped or aln.mapping_quality < args.min_mapq):
                continue
            qpos = pr.query_position
            if qpos is None or pr.is_refskip:
                continue
            if kind == "snv":
                if aln.query_sequence[qpos] != alt:
                    continue
                quals = aln.query_qualities
                if quals is not None and quals[qpos] < args.min_baseq:
                    continue
            elif kind == "del":
                if pr.indel != -(len(ref) - len(alt)):
                    continue
            else:  # ins
                n = len(alt) - 1
                if pr.indel != n or aln.query_sequence[qpos + 1:qpos + 1 + n] != alt[1:]:
                    continue
            hap = read_hap(aln)
            if hap is None:
                untagged += 1
            else:
                votes[hap] = votes.get(hap, 0) + 1
    return votes, untagged


def contig_votes(job):
    """Alt votes for one contig's records: [(record index, votes, untagged)]."""
    chrom, sites, args = job
    with pysam.AlignmentFile(args.alignments, reference_filename=args.reference) as bam:
        return [(idx, *alt_votes(bam, chrom, pos, ref, alt, kind, args)) for idx, pos, ref, alt, kind in sites]


def main(argv=None):
    args = parse_args(argv)
    with pysam.VariantFile(args.vcf) as vin:
        samples = list(vin.header.samples)
        if args.sample:
            if args.sample not in samples:
                sys.exit(f"phase_somatic_from_reads: sample {args.sample} not in {args.vcf} ({samples})")
            sample = args.sample
        elif len(samples) == 1:
            sample = samples[0]
        else:
            sys.exit(f"phase_somatic_from_reads: {args.vcf} has {len(samples)} samples, pass --sample")
        if "AHP" in vin.header.formats:
            sys.exit(f"phase_somatic_from_reads: {args.vcf} already declares FORMAT/AHP")

        with pysam.AlignmentFile(args.alignments, reference_filename=args.reference) as bam:
            contigs = set(bam.references)

        # First pass: the records to count reads for, grouped by contig
        sites = {}
        for idx, rec in enumerate(vin):
            if not has_alt(rec.samples[sample].get("GT")) or not rec.alts or len(rec.alts) != 1:
                continue
            kind = variant_kind(rec.ref, rec.alts[0])
            if kind is not None and rec.chrom in contigs:
                sites.setdefault(rec.chrom, []).append((idx, rec.pos, rec.ref, rec.alts[0], kind))

    jobs = [(chrom, recs, args) for chrom, recs in sites.items()]
    results = {}
    if args.threads > 1 and len(jobs) > 1:
        with multiprocessing.Pool(min(args.threads, len(jobs))) as pool:
            for part in pool.imap_unordered(contig_votes, jobs):
                results.update((idx, (votes, untagged)) for idx, votes, untagged in part)
    else:
        for job in jobs:
            results.update((idx, (votes, untagged)) for idx, votes, untagged in contig_votes(job))

    # Second pass: write every record, phased where the reads decide it
    vin = pysam.VariantFile(args.vcf)
    header = vin.header.copy()
    if "PS" not in header.formats:
        header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set: the germline phase block of the alt-supporting reads">')
    header.add_line('##FORMAT=<ID=AHP,Number=3,Type=Integer,Description="Alt-supporting reads on haplotype 1, on haplotype 2, and without a haplotag, in the haplotagged tumour alignments">')
    header.add_line(f"##phase_somatic_from_readsCommand={' '.join(sys.argv[1:] if argv is None else argv)}")
    vout = pysam.VariantFile(args.output, "wz", header=header)

    stats = dict(records=0, no_alt=0, unsupported=0, evaluated=0, not_het=0, too_few=0, ambiguous=0, phased=0)
    for idx, rec in enumerate(vin):
        stats["records"] += 1
        rec.translate(vout.header)
        call = rec.samples[sample]
        gt = call.get("GT")
        # Any phase from an earlier step is replaced by the read-based one
        if "PS" in call and call["PS"] is not None:
            call["PS"] = None
        if gt is not None:
            call.phased = False
            if is_het(gt):
                call["GT"] = (0, 1)
        if idx not in results:
            stats["unsupported" if has_alt(gt) else "no_alt"] += 1
            vout.write(rec)
            continue

        stats["evaluated"] += 1
        votes, untagged = results[idx]
        hp1 = sum(n for (_ps, hp), n in votes.items() if hp == 1)
        hp2 = sum(n for (_ps, hp), n in votes.items() if hp == 2)
        call["AHP"] = (hp1, hp2, untagged)
        if not is_het(gt):
            stats["not_het"] += 1
        elif hp1 + hp2 < args.min_reads:
            stats["too_few"] += 1
        else:
            (ps, hp), n = max(votes.items(), key=lambda kv: kv[1])
            if n / (hp1 + hp2) < args.min_fraction:
                stats["ambiguous"] += 1
            else:
                call["GT"] = (1, 0) if hp == 1 else (0, 1)
                call.phased = True
                call["PS"] = ps
                stats["phased"] += 1
        vout.write(rec)

    vout.close()
    vin.close()
    pysam.tabix_index(args.output, preset="vcf", force=True)
    print("\t".join(f"{k}={v}" for k, v in stats.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
