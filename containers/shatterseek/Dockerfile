# syntax=docker/dockerfile:1.6
FROM docker.io/mambaorg/micromamba:2.3.2

# R runtime for the ShatterSeek wrapper, with upstream parklab/ShatterSeek (not on conda) installed from a pinned commit
# plus a build-time patch adding T2T-CHM13 (patch_chm13.R). The wrapper ships with the pipeline in assets/shatterseek/.
COPY --chown=$MAMBA_USER containers/shatterseek/environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER containers/shatterseek/patch_chm13.R /tmp/patch_chm13.R

ARG SHATTERSEEK_COMMIT=3b3f4cddb23ff06a5d14d2ab5633d184e4f8d122
ADD --checksum=sha256:85d17017e144784f2bf8714c1746b91d2238d96bc3e5ae79195ae4c6bbeb7e4d --chown=$MAMBA_USER \
    https://github.com/parklab/ShatterSeek/archive/${SHATTERSEEK_COMMIT}.tar.gz /tmp/shatterseek.tar.gz
# T2T-CHM13 v2.0 cytobands, for the chm13 ideogram and chromosome sizes
ADD --checksum=sha256:09c7fa3ca3e222ebaabbecedc3086aa2543999824c0966f96b0c3f31c1126ea2 --chown=$MAMBA_USER \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed /tmp/chm13_cytobands.bed

RUN micromamba create --yes --name shatterseek --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && export PATH=/opt/conda/envs/shatterseek/bin:${PATH} \
    && mkdir -p /tmp/shatterseek-src \
    && tar -xzf /tmp/shatterseek.tar.gz -C /tmp/shatterseek-src --strip-components 1 \
    && (cd /tmp/shatterseek-src && Rscript /tmp/patch_chm13.R /tmp/chm13_cytobands.bed) \
    && R CMD INSTALL --no-docs --no-html /tmp/shatterseek-src \
    # smoke test: the bundled DO17373 example runs under every supported genome
    && Rscript -e 'suppressPackageStartupMessages(library(ShatterSeek)); data(DO17373); \
        for (g in c("hg19", "hg38", "chm13")) { set.seed(1); \
          r <- shatterseek(SV.sample = SVs(chrom1 = as.character(SV_DO17373$chrom1), pos1 = as.numeric(SV_DO17373$start1), \
                                           chrom2 = as.character(SV_DO17373$chrom2), pos2 = as.numeric(SV_DO17373$end2), \
                                           SVtype = as.character(SV_DO17373$svclass), strand1 = as.character(SV_DO17373$strand1), \
                                           strand2 = as.character(SV_DO17373$strand2)), \
                           seg.sample = CNVsegs(chrom = as.character(SCNA_DO17373$chromosome), start = SCNA_DO17373$start, \
                                                end = SCNA_DO17373$end, total_cn = SCNA_DO17373$total_cn), genome = g); \
          stopifnot(nrow(r@chromSummary) == 23); invisible(plot_chromothripsis(r, chr = "2", genome = g)) }; \
        cat("ShatterSeek", as.character(packageVersion("ShatterSeek")), "OK\n")' \
    && rm -rf /tmp/shatterseek.tar.gz /tmp/shatterseek-src /tmp/chm13_cytobands.bed /tmp/patch_chm13.R

ENV PATH=/opt/conda/envs/shatterseek/bin:${PATH} \
    SHATTERSEEK_COMMIT=${SHATTERSEEK_COMMIT}
