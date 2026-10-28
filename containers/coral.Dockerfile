# CoRAL with T2T-CHM13 support, built from a pinned commit of the fork.
#
# micromamba rather than python:slim because SCIP -- the open-source global MINLP
# solver CoRAL needs for its non-convex MIQCP -- is only packaged on conda-forge.
# Everything is installed into the base env and PATH is set explicitly, so the
# image works without shell activation (Nextflow/Singularity run no entrypoint).
FROM docker.io/mambaorg/micromamba:2.0.5-debian12-slim

ARG CORAL_REPO=https://github.com/robert-a-forsyth/CoRAL.git
ARG CORAL_REF=847f3d4
ARG CORAL_VERSION=3.0.0

LABEL org.opencontainers.image.title="CoRAL" \
      org.opencontainers.image.description="CoRAL amplicon reconstruction, T2T-CHM13 fork" \
      org.opencontainers.image.source="${CORAL_REPO}" \
      org.opencontainers.image.revision="${CORAL_REF}" \
      org.opencontainers.image.version="${CORAL_VERSION}" \
      org.opencontainers.image.licenses="BSD-3-Clause"

USER root
ENV PATH=/opt/conda/bin:$PATH \
    LANG=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1

# Build toolchain and the headers pysam and cvxopt fail without (upstream README)
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ git make pkg-config curl procps ca-certificates \
        libhdf5-dev libbz2-dev liblzma-dev zlib1g-dev \
        libcurl4-openssl-dev libssl-dev libffi-dev \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

# scip provides bin/scip, whose built-in AMPL reader is what Pyomo's SCIPAMPL
# plugin shells out to. pyscipopt is deliberately not installed -- it is unused.
# htslib is not needed: pysam installs from a manylinux wheel.
RUN micromamba install -y -n base -c conda-forge python=3.12 scip=10.1.0 \
    && micromamba clean --all --yes

# CPU-only torch first, or pomegranate/cnvkit pull the multi-GB CUDA wheel
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

RUN git clone "${CORAL_REPO}" /opt/CoRAL \
    && git -C /opt/CoRAL checkout "${CORAL_REF}" \
    && pip install --no-cache-dir /opt/CoRAL \
    && rm -rf /opt/CoRAL/.git

# Mount point only. Gurobi licences are user-supplied and never baked in.
RUN mkdir -p /opt/gurobi
ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

# No `| head`: SCIP dies on SIGPIPE when the reader closes early (exit 141)
RUN coral --help > /dev/null && scip -v > /dev/null

CMD ["coral", "--help"]
