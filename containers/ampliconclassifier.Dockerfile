# AmpliconClassifier v2.0.0 with T2T-CHM13 support, from a pinned fork commit.
#
# v2.0.0 is the first release with official CoRAL support. Installed with
# `pip install .` rather than a PATH symlink so the pinned BFBArchitect
# dependency, the console entry points and the bundled CHM13 lncRNA GFF3
# resource all land correctly.
FROM docker.io/library/python:3.11-slim

ARG AC_REPO=https://github.com/robert-a-forsyth/AmpliconClassifier.git
ARG AC_REF=cdeaa63
ARG AC_VERSION=2.0.0

LABEL org.opencontainers.image.title="AmpliconClassifier" \
      org.opencontainers.image.description="AmpliconClassifier 2.0.0, T2T-CHM13 fork" \
      org.opencontainers.image.source="${AC_REPO}" \
      org.opencontainers.image.revision="${AC_REF}" \
      org.opencontainers.image.version="${AC_VERSION}" \
      org.opencontainers.image.licenses="BSD-2-Clause"

ENV LANG=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ git procps ca-certificates \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        coinor-cbc \
    && rm -rf /var/lib/apt/lists/*

# BFBArchitect==1.0.1 comes in as a pinned dependency and pulls PuLP, CNVkit,
# pysam and matplotlib. gurobipy arrives too but needs no licence: BFBArchitect
# falls back Gurobi -> MOSEK -> CBC, so the image runs licence-free.
RUN pip install --no-cache-dir --upgrade pip \
    && git clone "${AC_REPO}" /opt/AmpliconClassifier \
    && git -C /opt/AmpliconClassifier checkout "${AC_REF}" \
    && pip install --no-cache-dir /opt/AmpliconClassifier \
    && rm -rf /opt/AmpliconClassifier/.git

# Mount point only. The AA data repo is ~1.1 GB of third-party-derived
# annotation with no stated licence, so it is staged at runtime.
RUN mkdir -p /opt/data_repo
ENV AA_DATA_REPO=/opt/data_repo

RUN amplicon_classifier.py --version \
    && python -c "import ampclasslib.ac_util as u; \
p = u.get_ncrna_file_loc('CHM13'); \
import os; assert os.path.exists(p), p; print('CHM13 lncRNA resource OK')"

CMD ["amplicon_classifier.py", "--help"]
