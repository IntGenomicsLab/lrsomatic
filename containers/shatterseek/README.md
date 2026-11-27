# ShatterSeek image

R runtime for the ShatterSeek wrapper used by the `SHATTERSEEK` module, with upstream
[parklab/ShatterSeek](https://github.com/parklab/ShatterSeek) (v1.1) installed from a pinned, checksum-verified commit
(`SHATTERSEEK_COMMIT` in the Dockerfile). Upstream supports hg19 and hg38 only; `patch_chm13.R` adds T2T-CHM13
(`genome = "chm13"`) at build time: the CHM13 v2.0 cytobands for the ideogram, the chromosome lengths for the
breakpoint-enrichment test, and `chm13` in the three genome checks. The statistics are upstream's, unchanged. The build
fails if the patch does not apply exactly, and runs the bundled DO17373 example under hg19, hg38 and chm13. The wrapper
is part of the pipeline in `assets/shatterseek/`. The module does not support `-profile conda`: the image ships the package.

Licence: ShatterSeek's `DESCRIPTION` declares GPL (>= 3), while the upstream README states that it is free for academic
use only.

Build and publish from the pipeline root (on the VSC login nodes: podman with the cgroupfs/file flags, TMPDIR on local
scratch). For a new upstream commit, update `SHATTERSEEK_COMMIT` and the `ADD --checksum` (the sha256 of
`https://github.com/parklab/ShatterSeek/archive/<commit>.tar.gz`) first. The SIF twin is what Singularity/Apptainer
pull, as for the ReConPlot and Padfoot images:

```bash
TAG=1.1-3b3f4cd-chm13
podman --cgroup-manager=cgroupfs --events-backend=file build --format docker \
    -f containers/shatterseek/Dockerfile -t "docker.io/robertaforsyth/shatterseek:$TAG" .
podman --cgroup-manager=cgroupfs --events-backend=file push "docker.io/robertaforsyth/shatterseek:$TAG"
# the SIF twin, built from the image just pushed and published with ORAS
# (`singularity remote login --username <user> oras://docker.io` first)
singularity build "shatterseek-$TAG.sif" "docker://docker.io/robertaforsyth/shatterseek:$TAG"
singularity push "shatterseek-$TAG.sif" "oras://docker.io/robertaforsyth/shatterseek-sif:$TAG"
```

The module pins both tags in `modules/local/shatterseek/main.nf`; override per site via
`process { withName: '.*:SHATTERSEEK_(SEVERUS_ASCAT|SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.
