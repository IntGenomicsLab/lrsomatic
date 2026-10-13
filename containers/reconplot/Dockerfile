# syntax=docker/dockerfile:1.6
FROM mambaorg/micromamba:2.3.2

# R runtime for the ReConPlot wrapper, with the ReConPlot package (not on conda) pre-installed
# from a pinned upstream commit. The wrapper itself ships with the pipeline in assets/reconplot/.
COPY --chown=$MAMBA_USER containers/reconplot/environment.yml /tmp/environment.yml

# Tarball checksum-verified like the Padfoot recipe (ADD --checksum needs the 1.6 Dockerfile frontend, Docker >= 25).
ARG RECONPLOT_COMMIT=88bb76011703d7908aaa2e5bdecfe18ba5737db8
ADD --checksum=sha256:9b5d19cb9a1ba5ce2b312549c076385aef15d3506a2085162c586ea179b018de --chown=$MAMBA_USER \
    https://github.com/cortes-ciriano-lab/ReConPlot/archive/${RECONPLOT_COMMIT}.tar.gz /tmp/reconplot.tar.gz

RUN micromamba create --yes --name reconplot --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && export PATH=/opt/conda/envs/reconplot/bin:${PATH} \
    && mkdir -p /tmp/reconplot-src \
    && tar -xzf /tmp/reconplot.tar.gz -C /tmp/reconplot-src --strip-components 1 \
    # Upstream bug: the per-chromosome clipping compares against the whole chr_selection vector instead of
    # the current row, which truncates chr1/chr2 in multi-chromosome (genome-wide) figures. Patch the two
    # lines to `chr_selection$chr[i]` and fail the build if the pattern is not found exactly twice.
    && test "$(grep -cF 'karyotype_data_now$chr == chr_selection$chr])' /tmp/reconplot-src/R/ReConPlot.R)" -eq 2 \
    && sed -i 's/karyotype_data_now\$chr == chr_selection\$chr\])/karyotype_data_now$chr == chr_selection$chr[i]])/g' /tmp/reconplot-src/R/ReConPlot.R \
    && test "$(grep -cF 'chr_selection$chr[i]])' /tmp/reconplot-src/R/ReConPlot.R)" -eq 2 \
    && R CMD INSTALL --no-docs --no-html /tmp/reconplot-src \
    && Rscript -e 'suppressPackageStartupMessages(library(ReConPlot)); cat("ReConPlot", as.character(packageVersion("ReConPlot")), "\n")' \
    && rm -rf /tmp/reconplot.tar.gz /tmp/reconplot-src

ENV PATH=/opt/conda/envs/reconplot/bin:${PATH}
