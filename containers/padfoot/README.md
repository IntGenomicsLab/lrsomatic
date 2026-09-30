# Padfoot RepeatMasker image

This is the image the `PADFOOT` module runs in (pinned by digest in `modules/local/padfoot/main.nf`). It contains:

- Padfoot itself at `/opt/padfoot` (`padfoot.py`, `padfoot/`, `beds/` with the hg38 and mm10 annotations), installed from
  the pinned commit of [Tim-Yu/Padfoot](https://github.com/Tim-Yu/Padfoot) given by `PADFOOT_COMMIT` in the Dockerfile,
  checksum-verified (Padfoot is not on bioconda; the fork adds SAVANA input support);
- its runtime dependencies from `containers/padfoot/environment.yml` (Python 3.12, pysam, pandas, biopython, samtools,
  minimap2, bedtools) and RepeatMasker 4.2.4;
- the Dfam 4.0 root and curated-consensus FamDB partitions (`FAMDB_DATA_DIR=/home/mambauser/dfam`), checksum-verified
  when the image is built.

The module does not support `-profile conda`: the image ships the tool, not just its dependencies.

Build and publish from the pipeline root:

```bash
export PADFOOT_IMAGE=<registry>/padfoot-repeatmasker:4.2.4-dfam4-padfoot-<padfoot commit>
docker build -f containers/padfoot/Dockerfile -t "$PADFOOT_IMAGE" .
docker push "$PADFOOT_IMAGE"
docker images --digests "$PADFOOT_IMAGE"   # the DIGEST column is what the module pins
```

Pin the pushed digest in the `container` directive of `modules/local/padfoot/main.nf`, or override per site:

```groovy
process { withName: '.*:PADFOOT_(SEVERUS_WAKHAN|SAVANA)' { container = '<registry>/padfoot-repeatmasker@sha256:<digest>' } }
```

To update Padfoot, change `PADFOOT_COMMIT` and the tarball checksum in the Dockerfile, rebuild, push and re-pin the digest.
The Dockerfile verifies the Padfoot tarball and the decompressed Dfam partition checksums, so a changed upstream file fails
the build rather than silently changing the image.
