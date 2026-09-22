# ReConPlot image

R runtime for the [ReConPlot wrapper](https://github.com/Tim-Yu/ReConPlot) used by the
`RECONPLOT` module, with the upstream [ReConPlot](https://github.com/cortes-ciriano-lab/ReConPlot)
R package (not distributed on conda) installed from a pinned commit. The wrapper scripts live in
`assets/reconplot/` and are staged by the pipeline.

Build and publish from the pipeline root:

```bash
export RECONPLOT_IMAGE=<registry>/reconplot:0.2-r4.4
docker build -f containers/reconplot/Dockerfile -t "$RECONPLOT_IMAGE" .
docker push "$RECONPLOT_IMAGE"
```

Pin the pushed digest in the `container` directive of `modules/local/reconplot/main.nf` (or override per site via
`process { withName: '.*:RECONPLOT_(SEVERUS_ASCAT|SEVERUS_WAKHAN|SAVANA)' { container = ... } }`). The module currently pins
`ghcr.io/tim-yu/reconplot@sha256:1145fc5aebe0227bec371f4c59b08b9a09871498e403c01b83f83973149ae9e7`.

Under `-profile conda` the module builds `modules/local/reconplot/environment.yml` and installs
ReConPlot at run time from the source tree staged via `--reconplot_pkg_url` / `--reconplot_pkg_dir`.
