# syntax=docker/dockerfile:1.6
FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER containers/padfoot/environment.yml /tmp/environment.yml

# Versioned Dfam directory: `releases/current` moves on with the next Dfam release and would break the rebuild
ARG DFAM_URL_BASE=https://www.dfam.org/releases/Dfam_4.0/families/FamDB
ARG DFAM_ROOT_SHA256=4cdb3d571126de3061d61c9e79d0ffcf0b1d9dc1d3c555bf9bcb8db695efa91a
ARG DFAM_CONSENSUS_SHA256=ce50a75fe90a5481864cd87792a80564667a4e7f373322730a42a7b737b8bbdd

RUN micromamba create --yes --name padfoot --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && export PATH=/opt/conda/envs/padfoot/bin:${PATH} \
    && mkdir -p /home/mambauser/dfam /tmp/repeatmasker-smoke \
    && wget --quiet --tries=3 --output-document=/tmp/dfam40.0.h5.gz "${DFAM_URL_BASE}/dfam40.0.h5.gz" \
    && wget --quiet --tries=3 --output-document=/tmp/dfam40.curated.consensus.0.h5.gz "${DFAM_URL_BASE}/dfam40.curated.consensus.0.h5.gz" \
    && gunzip /tmp/dfam40.0.h5.gz /tmp/dfam40.curated.consensus.0.h5.gz \
    && echo "${DFAM_ROOT_SHA256}  /tmp/dfam40.0.h5" | sha256sum --check \
    && echo "${DFAM_CONSENSUS_SHA256}  /tmp/dfam40.curated.consensus.0.h5" | sha256sum --check \
    && mv /tmp/dfam40.0.h5 /tmp/dfam40.curated.consensus.0.h5 /home/mambauser/dfam/ \
    && FAMDB_CONFIG="$(find /opt/conda/envs/padfoot/share -path '*/famdb-*/famdb.conf' -print -quit)" \
    && test -n "${FAMDB_CONFIG}" \
    && printf '[famdb]\nFAMDB_DATA_DIR = /home/mambauser/dfam\n' > "${FAMDB_CONFIG}" \
    && export PATH=/opt/conda/envs/padfoot/bin:${PATH} FAMDB_DATA_DIR=/home/mambauser/dfam \
    && famdb.py info \
    && printf '>smoke\nACGTACGTACGTACGT\n' > /tmp/repeatmasker-smoke/input.fa \
    && RepeatMasker -species human -dir /tmp/repeatmasker-smoke /tmp/repeatmasker-smoke/input.fa \
    && test -s /tmp/repeatmasker-smoke/input.fa.out \
    && rm -f /tmp/repeatmasker-smoke/input.fa.* \
    && RepeatMasker -species mus_musculus -dir /tmp/repeatmasker-smoke /tmp/repeatmasker-smoke/input.fa \
    && test -s /tmp/repeatmasker-smoke/input.fa.out \
    && rm -rf /tmp/dfam40.* /tmp/repeatmasker-smoke

# Padfoot itself is not on bioconda: install the pinned Tim-Yu/Padfoot commit (SAVANA input support) as source
# under /opt/padfoot, so the PADFOOT module needs no run-time download. Bump = new commit + checksum here, rebuild, push both images, bump the tags.
ARG PADFOOT_COMMIT=3846215096dc58f1e5c67084b74b56c857097885
ADD --checksum=sha256:61b9b10f70fd401562b57245480bb7fdf1dec82eb24e8bb15870e00bf8ca8c07 \
    https://github.com/Tim-Yu/Padfoot/archive/${PADFOOT_COMMIT}.tar.gz /tmp/padfoot.tar.gz

USER root
RUN mkdir -p /opt/padfoot \
    && tar -xzf /tmp/padfoot.tar.gz -C /opt/padfoot --strip-components 1 \
    && rm -f /tmp/padfoot.tar.gz \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/padfoot \
    && /opt/conda/envs/padfoot/bin/python3 /opt/padfoot/padfoot.py --version
USER $MAMBA_USER

ENV PATH=/opt/conda/envs/padfoot/bin:${PATH}
ENV FAMDB_DATA_DIR=/home/mambauser/dfam
ENV PADFOOT_HOME=/opt/padfoot
