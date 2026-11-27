# ReConPlot image

R runtime for the ReConPlot wrapper used by the `RECONPLOT` module, with the upstream
[ReConPlot](https://github.com/cortes-ciriano-lab/ReConPlot) R package (not distributed on conda)
installed from a pinned, checksum-verified commit (`RECONPLOT_COMMIT` in the Dockerfile). The Dockerfile
patches one upstream bug before installing: the per-chromosome clipping in `ReConPlot()` compared against the
whole `chr_selection$chr` vector instead of the current row, which truncated chr1/chr2 in multi-chromosome
(genome-wide) figures; the build fails if the two lines are not found. Images live on Docker Hub, like the
pipeline's other large custom images (`docker.io/timmy9527/reconplot`, SIF twin `docker.io/timmy9527/reconplot-sif`). The wrapper scripts are part of the
pipeline in `assets/reconplot/`. The module does not support `-profile conda`: the image ships the package.

Build and publish from the pipeline root:

```bash
TAG=0.2-r4.4-clipfix
docker build -f containers/reconplot/Dockerfile -t "docker.io/timmy9527/reconplot:$TAG" .
docker push "docker.io/timmy9527/reconplot:$TAG"
# the SIF twin, built from the image just pushed (by digest, so it is exactly the published image) and published with ORAS
# (`singularity remote login --username <user> oras://docker.io` first)
singularity build "reconplot-$TAG.sif" "docker://docker.io/timmy9527/reconplot:$TAG"
singularity push "reconplot-$TAG.sif" "oras://docker.io/timmy9527/reconplot-sif:$TAG"
```

The module pins both to the same tag in `modules/local/reconplot/main.nf` (`oras://docker.io/timmy9527/reconplot-sif:<tag>`
under Singularity/Apptainer, `docker.io/timmy9527/reconplot:<tag>` otherwise); override per site via
`process { withName: '.*:RECONPLOT_(SEVERUS_ASCAT|SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.

The R dependencies come from `containers/reconplot/environment.yml`. To update ReConPlot, change
`RECONPLOT_COMMIT` in the Dockerfile, rebuild, push both images and bump the tags.
