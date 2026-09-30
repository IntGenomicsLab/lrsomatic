# ReConPlot image

R runtime for the ReConPlot wrapper used by the `RECONPLOT` module, with the upstream
[ReConPlot](https://github.com/cortes-ciriano-lab/ReConPlot) R package (not distributed on conda)
installed from a pinned commit (`RECONPLOT_COMMIT` in the Dockerfile). The wrapper scripts are part of the
pipeline in `assets/reconplot/`. The module does not support `-profile conda`: the image ships the package.

Build and publish from the pipeline root:

```bash
TAG=0.2-r4.4
docker build -f containers/reconplot/Dockerfile -t "ghcr.io/tim-yu/reconplot:$TAG" .
docker push "ghcr.io/tim-yu/reconplot:$TAG"
# the SIF twin, built from the image just pushed and published with ORAS (singularity remote login first)
singularity build "reconplot-$TAG.sif" "docker-daemon://ghcr.io/tim-yu/reconplot:$TAG"
singularity push "reconplot-$TAG.sif" "oras://ghcr.io/tim-yu/reconplot-sif:$TAG"
```

The module pins both to the same tag in `modules/local/reconplot/main.nf` (`oras://ghcr.io/tim-yu/reconplot-sif:<tag>`
under Singularity/Apptainer, `ghcr.io/tim-yu/reconplot:<tag>` otherwise); override per site via
`process { withName: '.*:RECONPLOT_(SEVERUS_ASCAT|SEVERUS_WAKHAN|SAVANA)' { container = ... } }`.

The R dependencies come from `containers/reconplot/environment.yml`. To update ReConPlot, change
`RECONPLOT_COMMIT` in the Dockerfile, rebuild, push both images and bump the tags.
