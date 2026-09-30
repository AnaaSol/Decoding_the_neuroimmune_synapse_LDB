# Reproducible container for the Phase 1/1b/1c critical path (MAGMA gene-based
# test -> GCTA-COJO conditional analysis -> coloc.abf GWAS-vs-eQTL
# colocalization). Addresses the "renv pins versions but doesn't guarantee a
# clean OS-level rebuild" gap: this Dockerfile builds a self-contained image
# from a public base image with no dependency on this development machine's
# existing conda installation.
#
# Build:  docker build -t lbd-neuroimmune-phase1 .
# Run:    docker run --rm -v $(pwd)/data:/pipeline/data lbd-neuroimmune-phase1
#
# Scope: this image covers Phase 1/1b/1c (Python + R + MAGMA + GCTA + PLINK +
# tabix), the most genetics-heavy and reproducibility-sensitive part of the
# pipeline. Phases 2-4 depend on Seurat/BPCells/Harmony, which substantially
# increase image size and build time; they are not included in this image
# (see README's Docker section for the rationale on scoping this narrower).
#
# Uses conda-forge/bioconda for GCTA/PLINK/R (same packages verified working
# on the development machine) instead of hand-fetched binary URLs, which
# proved fragile (an initial attempt to wget a specific GCTA release URL
# 404'd/returned an invalid zip).

FROM condaforge/miniforge3:24.9.2-0

LABEL description="Decoding the Neuroimmune Synapse in LBD -- Phase 1/1b/1c (genetics + colocalization)"

RUN conda install -y -c conda-forge -c bioconda \
    r-base=4.3 r-recommended r-yaml \
    gcta=1.94.1 plink=1.90b6.21 htslib \
    libuv gmp \
    python=3.10 pandas numpy pyyaml \
    && conda clean -afy

# coloc + MendelianRandomization are not on conda-forge/bioconda in a
# pre-built form as of this Dockerfile's writing -- install from CRAN inside
# the image so the exact same R code path is exercised as on the host.
# r-recommended (Matrix/MASS/survival/etc.) and libuv/gmp (system libs for
# the fs/arrangements dependency chain) must be present BEFORE this step --
# an earlier attempt without them cascaded into ~6 failed packages.
RUN R -e "install.packages(c('coloc','MendelianRandomization'), repos='https://cloud.r-project.org')"

# Fail the build loudly if either package didn't actually install, instead
# of silently shipping a broken image (install.packages() alone does not
# make `docker build` exit non-zero on a single package failure).
RUN R -e "stopifnot(requireNamespace('coloc', quietly=TRUE), requireNamespace('MendelianRandomization', quietly=TRUE)); cat('coloc + MendelianRandomization OK\n')"

WORKDIR /pipeline
COPY scripts/phase1_genetic_prior/ /pipeline/scripts/phase1_genetic_prior/
COPY config/pipeline_config.yaml /pipeline/config/pipeline_config.yaml
COPY tools/magma /pipeline/tools/magma
RUN chmod +x /pipeline/tools/magma

# KNOWN LIMITATION (see README "Limitaciones conocidas"): every phase1
# script hardcodes project_dir as the absolute path on the development
# machine (/home/ana/Desktop/...) instead of reading it from
# scripts/utils/load_config.R, which exists but isn't actually wired into
# any script yet. Rather than silently patch each script's hardcoded path
# (or claim this is portable when it isn't), this symlink makes that literal
# host path resolve inside the container too, so the SAME unmodified scripts
# run identically in both places -- an honest workaround, not a fix.
RUN mkdir -p /home/ana/Desktop && ln -s /pipeline /home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB

# data/ is NOT baked into the image (public/re-derivable, too large -- mount
# as a volume at runtime: -v $(pwd)/data:/pipeline/data)
VOLUME ["/pipeline/data"]

CMD ["bash", "-c", "echo 'Mount data/ and run e.g.: Rscript scripts/phase1_genetic_prior/09_run_coloc.R brain_acc'"]
