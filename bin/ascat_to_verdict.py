#!/usr/bin/env python3
"""Rewrite ASCAT's purityploidy/segments tables into the two files ClairS-TO's tag_germline_variant.py reads.

Segment contigs are respelled to match the VCF (chr1 vs 1), since Verdict matches on exact contig name.
If ASCAT found no solution (purity 0), no purity file is written so Verdict tags nothing.
"""
import argparse
import gzip
import sys


def vcf_contigs(path):
    opener = gzip.open if path.endswith('.gz') else open
    contigs = []
    with opener(path, 'rt') as fp:
        for line in fp:
            if not line.startswith('#'):
                break
            if line.startswith('##contig=<'):
                fields = dict(kv.split('=', 1) for kv in line[len('##contig=<'):].rstrip('>\n').split(',') if '=' in kv)
                if 'ID' in fields:
                    contigs.append(fields['ID'])
    return contigs


def spell_like_vcf(chrom, contigs):
    """ASCAT's '1' becomes 'chr1' when that is what the VCF calls it, and the reverse."""
    if not contigs or chrom in contigs:
        return chrom
    if chrom.startswith('chr') and chrom[3:] in contigs:
        return chrom[3:]
    if 'chr' + chrom in contigs:
        return 'chr' + chrom
    return chrom


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--purityploidy', required=True, help="ASCAT <sample>.purityploidy.txt")
    parser.add_argument('--segments', required=True, help="ASCAT <sample>.segments.txt")
    parser.add_argument('--vcf', required=True, help="a ClairS-TO VCF; only its ##contig header lines are read")
    parser.add_argument('--sample', required=True, help="value of the Sample column")
    parser.add_argument('--purity_out', required=True)
    parser.add_argument('--cna_out', required=True)
    args = parser.parse_args()

    with open(args.purityploidy) as fp:
        header = fp.readline().rstrip('\n').split('\t')
        values = fp.readline().rstrip('\n').split('\t')
    try:
        purity = float(values[header.index('AberrantCellFraction')])
        ploidy = float(values[header.index('Ploidy')])
    except (ValueError, IndexError):
        sys.exit(f"[ERROR] {args.purityploidy} is not an ASCAT purityploidy table: {header} / {values}")

    contigs = vcf_contigs(args.vcf)
    unmatched = set()
    n_segments = 0
    with open(args.segments) as fp, open(args.cna_out, 'w') as out:
        header = fp.readline().rstrip('\n').split('\t')
        col = {name: header.index(name) for name in ('chr', 'startpos', 'endpos', 'nMajor', 'nMinor')}
        out.write('Sample\tChromosome\tStartPosition\tEndPosition\tnMajor\tnMinor\n')
        for line in fp:
            fields = line.rstrip('\n').split('\t')
            chrom = spell_like_vcf(fields[col['chr']], contigs)
            if contigs and chrom not in contigs:
                unmatched.add(chrom)
            out.write('\t'.join([args.sample, chrom, fields[col['startpos']], fields[col['endpos']],
                                 fields[col['nMajor']], fields[col['nMinor']]]) + '\n')
            n_segments += 1
    if unmatched:
        print(f"[WARNING] {len(unmatched)} ASCAT contig(s) are not in the VCF header and will match no variant: "
              f"{', '.join(sorted(unmatched))}", file=sys.stderr)

    if not purity > 0:
        print(f"[WARNING] ASCAT found no purity/ploidy solution ({purity}, {ploidy}); "
              f"not writing {args.purity_out}, so Verdict tags nothing.", file=sys.stderr)
        return
    with open(args.purity_out, 'w') as out:
        out.write('Sample\tPurity\tPloidy\tGoodnessOfFit\n')
        out.write(f'{args.sample}\t{purity}\t{ploidy}\tNA\n')
    print(f"[INFO] purity {purity}, ploidy {ploidy}, {n_segments} segments written for {args.sample}")


if __name__ == '__main__':
    main()
