# Padfoot RepeatMasker image

This is the image the `PADFOOT` module runs in (pinned by digest in `modules/local/padfoot/main.nf`). It embeds the Dfam 4.0 root and curated-consensus FamDB partitions, configures FamDB at build time, and validates human RepeatMasker during the build.

Build and publish from the pipeline root:

```bash
export PADFOOT_IMAGE=<registry>/padfoot-repeatmasker:4.2.4-dfam4
docker build -f containers/padfoot/Dockerfile -t "$PADFOOT_IMAGE" .
docker push "$PADFOOT_IMAGE"
```

Resolve the pushed digest and pin it in the `container` directive of `modules/local/padfoot/main.nf`, or override per site:

```groovy
process { withName: '.*:PADFOOT_(SEVERUS_WAKHAN|SAVANA)' { container = '<registry>/padfoot-repeatmasker@sha256:<digest>' } }
```

The Dockerfile verifies the decompressed Dfam partition checksums. A changed Dfam `current` release therefore fails the build rather than silently changing the annotation database.

The module currently pins `ghcr.io/tim-yu/padfoot-repeatmasker@sha256:f98b0d352ec47cd9fa015f321959dcba7f6c6cd4da05beaea8d811b97dff70f5`; RepeatMasker runs by default (`--padfoot_run_repeatmasker`).
