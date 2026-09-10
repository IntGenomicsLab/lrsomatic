FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER modules/local/padfoot/environment.yml /tmp/environment.yml

ARG DFAM_URL_BASE=https://www.dfam.org/releases/current/families/FamDB
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
    && rm -rf /tmp/dfam40.* /tmp/repeatmasker-smoke

ENV PATH=/opt/conda/envs/padfoot/bin:${PATH}
ENV FAMDB_DATA_DIR=/home/mambauser/dfam
