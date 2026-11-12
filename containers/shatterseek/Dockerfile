# syntax=docker/dockerfile:1.6
FROM docker.io/mambaorg/micromamba:2.3.2

# R runtime for the ShatterSeek wrapper, with the robert-a-forsyth/ShatterSeek fork (chm13 support, statistical
# criteria fixed to the PCAWG definitions, PCAWG confidence tiers) installed from a pinned commit. The wrapper ships
# with the pipeline in assets/shatterseek/.
COPY --chown=$MAMBA_USER containers/shatterseek/environment.yml /tmp/environment.yml

ARG SHATTERSEEK_COMMIT=e1cee6c8ee75ea1a1e9f4d1c57b7a8a0a199ae91
ADD --checksum=sha256:74b4258a0e6b509aeba4976348feebf910c098b9556e85008795d1b4e6286b69 --chown=$MAMBA_USER \
    https://github.com/robert-a-forsyth/ShatterSeek/archive/${SHATTERSEEK_COMMIT}.tar.gz /tmp/shatterseek.tar.gz

RUN micromamba create --yes --name shatterseek --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && export PATH=/opt/conda/envs/shatterseek/bin:${PATH} \
    && mkdir -p /tmp/shatterseek-src \
    && tar -xzf /tmp/shatterseek.tar.gz -C /tmp/shatterseek-src --strip-components 1 \
    && R CMD INSTALL --no-docs --no-html /tmp/shatterseek-src \
    # the fork's regression test: DO17373 against the PCAWG-published values, plus synthetic hg38 / chm13 events
    && Rscript /tmp/shatterseek-src/tests/test_do17373.R \
    && rm -rf /tmp/shatterseek.tar.gz /tmp/shatterseek-src

ENV PATH=/opt/conda/envs/shatterseek/bin:${PATH} \
    SHATTERSEEK_COMMIT=${SHATTERSEEK_COMMIT}
