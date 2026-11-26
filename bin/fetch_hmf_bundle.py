#!/usr/bin/env python3
"""Download, MD5-check and unpack a hosted HMF report bundle (one tarball per genome).

The bundle is an oncoanalyser-layout directory (common/, dna/sv/, dna/variants/) packed as a single
top-level directory; it is unpacked without that top level into --outdir, so the report's slots are at
fixed paths (common/DriverGenePanel.38.tsv, common/ensembl_data/, ...). Standard library only.
"""

import argparse
import hashlib
import os
import shutil
import sys
import tarfile
import time
import urllib.request

REQUIRED = [
    "common/DriverGenePanel.38.tsv",
    "common/ensembl_data/ensembl_gene_data.csv",
    "common/ensembl_data/ensembl_trans_exon_data.csv",
]


def download(url, dest, tries=10, wait=10):
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as out:
                shutil.copyfileobj(r, out)
            return
        except Exception as e:  # network hiccups: retry, as wget --tries does elsewhere
            if attempt == tries:
                sys.exit(f"ERROR: could not download {url}: {e}")
            sys.stderr.write(f"download attempt {attempt} failed ({e}); retrying in {wait}s\n")
            time.sleep(wait)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", required=True)
    ap.add_argument("--md5", default="", help="Expected MD5 of the tarball; empty skips the check")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    archive = os.path.basename(a.url.split("?")[0])
    download(a.url, archive)
    if a.md5:
        h = hashlib.md5()
        with open(archive, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != a.md5:
            sys.exit(f"ERROR: {archive} has MD5 {h.hexdigest()}, expected {a.md5}")

    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        tops = {m.name.split("/", 1)[0] for m in members}
        if len(tops) != 1:
            sys.exit(f"ERROR: {archive} must hold one top-level directory, found {sorted(tops)}")
        top = tops.pop()
        for m in members:
            rel = m.name[len(top):].lstrip("/")
            if not rel or rel.startswith("..") or os.path.isabs(rel) or not (m.isfile() or m.isdir()):
                continue
            m.name = rel
            tar.extract(m, a.outdir, filter="data")
    os.remove(archive)

    missing = [p for p in REQUIRED if not os.path.exists(os.path.join(a.outdir, p))]
    if missing:
        sys.exit(f"ERROR: {archive} lacks {', '.join(missing)}")


if __name__ == "__main__":
    main()
