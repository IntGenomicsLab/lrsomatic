#!/usr/bin/env python3
"""Convert ASCAT's cnvs.txt into the headerless BED CoRAL's --cn-seg expects.

CoRAL reads total copy number from the last column and builds chromosome sizes from
the BAM header, so contigs are respelled to match the reference (chr1 vs 1).
"""
import argparse
import sys

REQUIRED = ('chr', 'startpos', 'endpos', 'nMajor', 'nMinor')


def fai_contigs(path):
    if not path:
        return []
    with open(path) as fp:
        return [line.split('\t', 1)[0] for line in fp if line.strip()]


def spell_like_reference(chrom, contigs):
    """ASCAT's '1' becomes 'chr1' when that is what the reference calls it, and the reverse."""
    if not contigs or chrom in contigs:
        return chrom
    if chrom.startswith('chr') and chrom[3:] in contigs:
        return chrom[3:]
    if 'chr' + chrom in contigs:
        return 'chr' + chrom
    return chrom


def sort_key(chrom):
    """Natural contig order: 1-22, then X, Y, then anything else alphabetically."""
    bare = chrom[3:] if chrom.startswith('chr') else chrom
    if bare.isdigit():
        return (0, int(bare), '')
    if bare in ('X', 'Y'):
        return (1, 'XY'.index(bare), '')
    return (2, 0, bare)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cnvs', required=True, help="ASCAT <sample>.cnvs.txt")
    parser.add_argument('--fai', help="Reference .fai, to respell contigs to match the BAM")
    parser.add_argument('--output', required=True, help="BED4 written for CoRAL --cn-seg")
    args = parser.parse_args()

    contigs = fai_contigs(args.fai)

    with open(args.cnvs) as fp:
        header = fp.readline().rstrip('\n').split('\t')
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            sys.exit(f"ERROR: {args.cnvs} is missing required columns: {', '.join(missing)}")
        idx = {c: header.index(c) for c in REQUIRED}

        rows, dropped = [], 0
        for line in fp:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            start, end = int(fields[idx['startpos']]), int(fields[idx['endpos']])
            # CoRAL's segment parser rejects these, so drop them here with a count
            if start >= end:
                dropped += 1
                continue
            total_cn = round(float(fields[idx['nMajor']])) + round(float(fields[idx['nMinor']]))
            rows.append((spell_like_reference(fields[idx['chr']], contigs), start, end, total_cn))

    if dropped:
        print(f"WARNING: dropped {dropped} segments with start >= end", file=sys.stderr)

    rows.sort(key=lambda r: (sort_key(r[0]), r[1]))
    with open(args.output, 'w') as out:
        for chrom, start, end, total_cn in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{total_cn}\n")

    if not rows:
        print(f"WARNING: {args.output} is empty; CoRAL will find no seeds", file=sys.stderr)


if __name__ == '__main__':
    main()
