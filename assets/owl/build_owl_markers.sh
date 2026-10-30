#!/usr/bin/env bash
# Build an owl marker BED from a RepeatMasker .out file, following owl's analysis/MARKERS.md recipe.
# Usage: build_owl_markers.sh <repeatmasker.out[.gz]> <output.bed>
set -euo pipefail

OWL_SHA=ad60d78603b6a477200ae753b28e1c5af985c35e
here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

curl -fsSL "https://raw.githubusercontent.com/PacificBiosciences/owl/${OWL_SHA}/scripts/repfilter.py" -o "$tmp/repfilter.py"
python3 "$here/rmsk_out_to_txt.py" "$1" > "$tmp/rmsk.txt"
python3 "$tmp/repfilter.py" --min-length 20 --max-length 100 --min-identity 97 --max-unit-len 6 "$tmp/rmsk.txt" \
    | grep -v -E "Un|random|alt|fix|chrY" \
    | cut -f 1-4 \
    | LC_ALL=C sort -k1,1 -k2,2n > "$2"
