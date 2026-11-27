#!/usr/bin/env python3
"""Convert ASCAT's cnvs.txt into the headerless BED CoRAL's --cn-seg expects.

ASCAT segments are 1-based inclusive; the BED is 0-based half-open and tiles each
chromosome without gaps (gaps split at their midpoint, ends run to the contig edges).
CoRAL reads total copy number from the last column and builds chromosome sizes from
the BAM header, so contigs are respelled to match the reference (chr1 vs 1).
"""
import argparse
import sys

REQUIRED = ('chr', 'startpos', 'endpos', 'nMajor', 'nMinor')


def fai_lengths(path):
    """Contig name -> length, from a .fai; empty without one."""
    if not path:
        return {}
    lengths = {}
    with open(path) as fp:
        for line in fp:
            if line.strip():
                fields = line.split('\t')
                lengths[fields[0]] = int(fields[1])
    return lengths


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


def tile(segments, length):
    """Make one chromosome's sorted [start, end, cn] segments contiguous from 0 to its length.

    Returns the tiles and how many segments lay wholly inside an earlier one and were dropped.
    """
    tiles, contained = [], 0
    for start, end, cn in segments:
        if not tiles:
            tiles.append([0, end, cn])
            continue
        prev = tiles[-1]
        if end <= prev[1]:
            contained += 1
            continue
        # The midpoint of a gap, or of an overlap; both neighbours keep at least 1 bp
        boundary = min(max((prev[1] + start) // 2, prev[0] + 1), end - 1)
        prev[1] = boundary
        tiles.append([boundary, end, cn])
    if tiles and length and length > tiles[-1][0]:
        tiles[-1][1] = length
    return tiles, contained


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cnvs', required=True, help="ASCAT <sample>.cnvs.txt")
    parser.add_argument('--fai', help="Reference .fai, to respell contigs and extend the last segment to each contig's end")
    parser.add_argument('--output', required=True, help="BED4 written for CoRAL --cn-seg")
    args = parser.parse_args()

    lengths = fai_lengths(args.fai)

    with open(args.cnvs) as fp:
        header = fp.readline().rstrip('\n').split('\t')
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            sys.exit(f"ERROR: {args.cnvs} is missing required columns: {', '.join(missing)}")
        idx = {c: header.index(c) for c in REQUIRED}

        by_chrom, dropped = {}, 0
        for line in fp:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            start, end = int(fields[idx['startpos']]), int(fields[idx['endpos']])
            # start == end is a valid one-base segment; only an inverted one is malformed
            if start > end:
                dropped += 1
                continue
            total_cn = round(float(fields[idx['nMajor']])) + round(float(fields[idx['nMinor']]))
            chrom = spell_like_reference(fields[idx['chr']], lengths)
            by_chrom.setdefault(chrom, []).append((start - 1, end, total_cn))

    if dropped:
        print(f"WARNING: dropped {dropped} segments with start > end", file=sys.stderr)

    contained = 0
    with open(args.output, 'w') as out:
        for chrom in sorted(by_chrom, key=sort_key):
            tiles, n = tile(sorted(by_chrom[chrom]), lengths.get(chrom))
            contained += n
            for start, end, total_cn in tiles:
                out.write(f"{chrom}\t{start}\t{end}\t{total_cn}\n")

    if contained:
        print(f"WARNING: dropped {contained} segments lying wholly inside another", file=sys.stderr)
    if not by_chrom:
        print(f"WARNING: {args.output} is empty; CoRAL will find no seeds", file=sys.stderr)


if __name__ == '__main__':
    main()
