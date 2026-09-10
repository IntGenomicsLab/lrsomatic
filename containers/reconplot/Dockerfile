FROM mambaorg/micromamba:2.3.2

# R runtime for the ReConPlot wrapper, with the ReConPlot package (not on conda) pre-installed
# from a pinned upstream commit. The wrapper itself is staged by the pipeline (params.reconplot_url).
COPY --chown=$MAMBA_USER modules/local/reconplot/environment.yml /tmp/environment.yml

ARG RECONPLOT_COMMIT=88bb76011703d7908aaa2e5bdecfe18ba5737db8
ADD --chown=$MAMBA_USER https://github.com/cortes-ciriano-lab/ReConPlot/archive/${RECONPLOT_COMMIT}.tar.gz /tmp/reconplot.tar.gz

RUN micromamba create --yes --name reconplot --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && export PATH=/opt/conda/envs/reconplot/bin:${PATH} \
    && mkdir -p /tmp/reconplot-src \
    && tar -xzf /tmp/reconplot.tar.gz -C /tmp/reconplot-src --strip-components 1 \
    && R CMD INSTALL --no-docs --no-html /tmp/reconplot-src \
    && Rscript -e 'suppressPackageStartupMessages(library(ReConPlot)); cat("ReConPlot", as.character(packageVersion("ReConPlot")), "\n")' \
    && rm -rf /tmp/reconplot.tar.gz /tmp/reconplot-src

ENV PATH=/opt/conda/envs/reconplot/bin:${PATH}
