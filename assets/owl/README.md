# owl marker sets

owl ([PacificBiosciences/owl](https://github.com/PacificBiosciences/owl)) scores microsatellite instability (MSI) at a fixed set of marker loci: a 4-column BED of chromosome, start, end and motif written as `(AC)n`.

- **GRCh38:** PacBio's curated `data/GRCh38_owl_markers.bed.gz`, fetched from the owl repository at commit `ad60d78` (see `conf/igenomes.config`).
- **CHM13:** no marker set is published (pbbioconda issue #810), so this one is built here with owl's own recipe (`analysis/MARKERS.md`) from UCSC's hs1 RepeatMasker annotation.

## Building the CHM13 set

```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz   # md5 5cc807e86912c46af52d3bd35f2b833c
assets/owl/build_owl_markers.sh hs1.repeatMasker.out.gz CHM13_owl_markers.bed
bgzip CHM13_owl_markers.bed && tabix -p bed CHM13_owl_markers.bed.gz
```

- **`rmsk_out_to_txt.py`** converts RepeatMasker's `.out` layout to the 17-column UCSC `rmsk.txt` layout that owl's `scripts/repfilter.py` reads.
- **`build_owl_markers.sh`** fetches `repfilter.py` at the pinned owl commit, then applies the recipe unchanged: repeats of 20-100 bp, at least 97% identity, motifs up to 6 bp, and only chr1-22 and chrX.

The result is 170,612 loci (`CHM13_owl_markers.bed.gz` md5 `d0259b2d5fef9df595741c3c8b71abf4`), hosted on the IntGenomicsLab Globus collection.

## Checks

- **Converter:** converting UCSC's hg38 `hg38.fa.out.gz` gives the same rows as UCSC's hg38 `rmsk.txt.gz` for 5,520,017 of 5,520,118 repeats (all differences are SINE/LINE/DNA rows). After the recipe the two give the identical 164,374 loci.
- **Curation:** PacBio's curated GRCh38 set (146,562 loci) is a strict subset of the recipe's hg38 output (164,374 loci). Its extra pruning step, for no-calls and population variability, removed 17,812 loci (10.8%), and is not public. The CHM13 set has no such pruning, so CHM13 `%high` is not on the same scale as GRCh38 `%high`.
