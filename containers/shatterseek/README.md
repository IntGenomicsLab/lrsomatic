# ShatterSeek image

R runtime for the ShatterSeek wrapper used by the `SHATTERSEEK` module, with the
[robert-a-forsyth/ShatterSeek](https://github.com/robert-a-forsyth/ShatterSeek/tree/chm13-support) fork installed from
a pinned, checksum-verified commit (`SHATTERSEEK_COMMIT` in the Dockerfile). Relative to upstream
[parklab/ShatterSeek](https://github.com/parklab/ShatterSeek) the fork adds T2T-CHM13 (`genome = "chm13"`), corrects
the statistical criteria to the PCAWG definitions, keeps one cluster per chromosome and labels clusters with the PCAWG
confidence tiers (`classify_chromothripsis()`); its README lists every change. The build runs the fork's regression
test (`tests/test_do17373.R`: the DO17373 example against the PCAWG-published values, plus synthetic hg38 and chm13
events) and fails if it does not pass. The wrapper is part of the pipeline in `assets/shatterseek/`. The module does
not support `-profile conda`: the image ships the package.

Licence: the fork's `DESCRIPTION` declares GPL (>= 3), as upstream does, while the upstream README states that
ShatterSeek is free for academic use only.

Build and publish from the pipeline root (on the VSC login nodes: podman with the cgroupfs/file flags, TMPDIR on local
scratch). After a new fork commit, update `SHATTERSEEK_COMMIT` and the `ADD --checksum` (the sha256 of
`https://github.com/robert-a-forsyth/ShatterSeek/archive/<commit>.tar.gz`) first:

```bash
TAG=1.1.1-chm13-<short commit>
podman --cgroup-manager=cgroupfs --events-backend=file build --format docker \
    -f containers/shatterseek/Dockerfile -t "docker.io/robertaforsyth/shatterseek:$TAG" .
podman --cgroup-manager=cgroupfs --events-backend=file push "docker.io/robertaforsyth/shatterseek:$TAG"
```

There is no SIF twin: like CoRAL, Singularity/Apptainer pull the Docker image. The module pins the tag in
`modules/local/shatterseek/main.nf`; override per site via
`process { withName: '.*:SHATTERSEEK_(SEVERUS_ASCAT|SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.
