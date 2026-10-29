#!/usr/bin/env python3
"""Drop SAVANA SVs that match a germline SV panel of normals (PoN).

SAVANA's PacBio classifier (`savana to --pb`) has no tumour-only mode: it rejects a call
only on normal-read evidence, so without a normal every germline SV with enough tumour
reads is written to the somatic VCF (~20-30k per genome). This filters that VCF against
the Severus PoN the pipeline already uses, with Severus's own matching rule (extract_pon
and add_pon in severus/breakpoint_finder.py), so both callers call germline the same way.

PoN rows: chr1,pos1,chr2,pos2,ci1,ci2,svtype,fraction -- for INS, pos2 is the length.
A header row (the CHM13 PoN has one) is skipped because its pos1 is not a number.

An event is matched on its lower breakpoint, and both BND mates are kept or dropped together.
"""
import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict

# Severus add_pon constants
BUFF = 2000
CLUST_LEN = 250
MAX_LEN_DIFF = 50
MAX_SUM_DIFF = 1000

MATE_RE = re.compile(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_pon(path):
    """Key (chrom, or (chrom1, chrom2) between chromosomes) -> sorted columns pos1, ci1, end, ci2, partner chrom."""
    rows = defaultdict(list)
    with open_text(path) as fp:
        for line in fp:
            fields = line.strip().split(',')
            if len(fields) < 7:
                continue
            c1, p1, c2, p2, ci1, ci2, svtype = fields[:7]
            try:
                p1, p2, ci1, ci2 = int(p1), int(p2), int(ci1), int(ci2)
            except ValueError:
                continue
            if svtype == 'INS':
                rows[c1].append((p1, ci1, p1 + p2, ci2, c1))
            elif c1 == c2:
                rows[c1].append((p1, ci1, p2, ci2, c2))
            else:
                rows[(c1, c2)].append((p1, ci1, p2, ci2, c2))
                rows[(c2, c1)].append((p2, ci2, p1, ci1, c1))
    pon = {}
    for key, entries in rows.items():
        entries.sort()
        pon[key] = list(zip(*entries))
    return pon


def info_dict(info):
    out = {}
    for kv in info.split(';'):
        key, _, value = kv.partition('=')
        out[key] = value
    return out


def breakpoints(fields):
    """chrom, pos, partner chrom, partner pos, length, is_ins for one VCF record."""
    chrom, pos, alt = fields[0], int(fields[1]), fields[4]
    info = info_dict(fields[7])
    if info.get('BP_NOTATION') == '<INS>' or info.get('SVTYPE') == 'INS' or alt == '<INS>':
        try:
            length = abs(int(float(info.get('SVLEN', '0'))))
        except ValueError:
            length = 0
        return chrom, pos, chrom, pos, length, True
    mate = MATE_RE.search(alt)
    if mate:
        chrom2, pos2 = mate.group(1), int(mate.group(2))
    else:
        chrom2, pos2 = chrom, int(info.get('END', pos))
    return chrom, pos, chrom2, pos2, abs(pos2 - pos), False


def in_pon(pon, chrom, pos, chrom2, pos2, length):
    """Severus add_pon, with no breakpoint confidence interval for the query (SAVANA reports none)."""
    entries = pon.get(chrom if chrom == chrom2 else (chrom, chrom2))
    if not entries:
        return False
    starts, cis1, ends, cis2, partners = entries
    lo = bisect.bisect_left(starts, pos - BUFF)
    hi = bisect.bisect_right(starts, pos + BUFF)
    for i in range(lo, hi):
        if partners[i] != chrom2:
            continue
        if chrom == chrom2:
            len_diff = abs(ends[i] - starts[i] - length)
        else:
            len_diff = abs(ends[i] - pos2)
        sum_diff = abs(starts[i] - pos) + len_diff
        if (len_diff <= MAX_LEN_DIFF and sum_diff <= MAX_SUM_DIFF) or sum_diff <= max(cis1[i], cis2[i]) + CLUST_LEN:
            return True
    return False


def event_id(record_id):
    return re.sub(r'_[12]$', '', record_id)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='SAVANA classified somatic VCF')
    parser.add_argument('--pon', required=True, help='Severus SV panel of normals (.tsv or .tsv.gz)')
    parser.add_argument('--output', required=True, help='Filtered VCF')
    parser.add_argument('--summary', required=True, help='TSV with the event counts before and after')
    parser.add_argument('--sample', default='', help='Sample name for the summary')
    args = parser.parse_args()

    pon = read_pon(args.pon)

    header, records = [], []
    with open_text(args.vcf) as fp:
        for line in fp:
            if line.startswith('#'):
                header.append(line)
            elif line.strip():
                records.append(line)

    # one decision per event, made on its lower breakpoint (Severus's bp_1)
    lowest = {}
    for line in records:
        fields = line.rstrip('\n').split('\t')
        bp = breakpoints(fields)
        eid = event_id(fields[2])
        if eid not in lowest or (bp[0], bp[1]) < (lowest[eid][0], lowest[eid][1]):
            lowest[eid] = bp
    germline = {eid for eid, (c, p, c2, p2, length, _ins) in lowest.items() if in_pon(pon, c, p, c2, p2, length)}

    kept = 0
    with open(args.output, 'w') as out:
        for line in header:
            if line.startswith('#CHROM'):
                out.write('##savana_pon_filter=<Description="Events matching the germline SV panel of normals '
                          f'{args.pon.split("/")[-1]} removed (Severus PoN matching rule)">\n')
            out.write(line)
        for line in records:
            if event_id(line.split('\t', 3)[2]) not in germline:
                out.write(line)
                kept += 1

    n_events = len(lowest)
    with open(args.summary, 'w') as fp:
        fp.write('sample\tevents_in\tevents_pon_matched\tevents_out\trecords_in\trecords_out\n')
        fp.write(f'{args.sample}\t{n_events}\t{len(germline)}\t{n_events - len(germline)}\t{len(records)}\t{kept}\n')
    print(f'{args.sample}: {len(germline)} of {n_events} SAVANA events match the PoN and were removed; '
          f'{n_events - len(germline)} kept ({kept} of {len(records)} records)', file=sys.stderr)


if __name__ == '__main__':
    main()
