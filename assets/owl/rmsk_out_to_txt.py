#!/usr/bin/env python3
"""Convert a RepeatMasker .out file to the 17-column UCSC rmsk.txt layout read by owl's repfilter.py."""

import gzip
import sys


def open_text(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def num(field):
    # rmsk.txt stores parenthesised .out counts as negatives
    return -int(field.strip("()")) if field.startswith("(") else int(field)


def convert(line):
    f = line.split()
    # Data rows have 15 fields, plus a trailing '*' when a better-scoring repeat overlaps
    if len(f) not in (15, 16) or not f[0].isdigit():
        return None
    score, div, dele, ins, chrom, qbeg, qend, qleft, strand, name, cls_fam = f[:11]
    cls, _, fam = cls_fam.partition("/")
    return [
        0, int(score),
        round(float(div) * 10), round(float(dele) * 10), round(float(ins) * 10),
        chrom, int(qbeg) - 1, int(qend), num(qleft),
        "-" if strand == "C" else "+",
        name, cls, fam or cls,
        num(f[11]), num(f[12]), num(f[13]), int(f[14]),
    ]


def main():
    with open_text(sys.argv[1] if len(sys.argv) > 1 else "-") as fh:
        for line in fh:
            row = convert(line)
            if row:
                print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
