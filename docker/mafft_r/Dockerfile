# Dockerfile: R + MAFFT container for use with Nextflow
#
# Base image: rocker/tidyverse gives us R + tidyverse pre-installed and
# pre-compiled, which saves a lot of build time vs. installing tidyverse
# from scratch (many of its dependencies need compilation).
FROM rocker/tidyverse:4.4.1

LABEL description="R + Bioconductor + MAFFT container for Nextflow pipelines"
LABEL maintainer="you@example.com"

# -----------------------------------------------------------------------
# System dependencies
# -----------------------------------------------------------------------
# - mafft: multiple sequence alignment tool
# - procps: provides `ps`, which Nextflow needs inside containers to poll
#           process status/resource usage
# - locales/build tools: common requirements for R package compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    mafft \
    procps \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------
# R packages
# -----------------------------------------------------------------------
# tidyverse is already included in the rocker/tidyverse base image.
# Biostrings comes from Bioconductor; phylotools and readxl come from CRAN.
#
# IMPORTANT: install.packages() does NOT make R exit with a non-zero status
# just because a package failed to compile/install -- it only emits a
# warning. That means a naive `RUN R -e "install.packages(...)"` can "pass"
# in Docker even though the package silently never installed. To avoid that,
# every install below is followed by an explicit requireNamespace() check
# that calls quit(status = 1) if the package isn't actually usable, so the
# build fails immediately and clearly instead of failing later downstream.
#
# Parallelism, two independent knobs:
#   - MAKEFLAGS="-j$(nproc)": lets each individual package's C/Fortran source
#     compile using all cores (matters for e.g. ape, which phylotools needs).
#     Must be set on the RUN line itself (or exported first) -- it can't be
#     set with a Dockerfile ENV instruction, since ENV doesn't expand
#     $(nproc); only a shell (which RUN uses) does that substitution.
#   - Ncpus = parallel::detectCores() inside install.packages(): lets R build
#     several DIFFERENT packages concurrently in one call. This is what
#     actually pays off for the phylotools call below, since
#     dependencies = TRUE pulls in ape and friends alongside it.
# Binary installs from a package manager (like Posit's, which CRAN already
# points to in rocker/tidyverse) skip compilation entirely, so neither knob
# does anything for those -- they're only relevant for source builds.

RUN MAKEFLAGS="-j$(nproc)" R -e "\
    options(Ncpus = parallel::detectCores()); \
    install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    if (!requireNamespace('BiocManager', quietly = TRUE)) quit(status = 1)"

RUN MAKEFLAGS="-j$(nproc)" R -e "\
    options(Ncpus = parallel::detectCores()); \
    BiocManager::install('Biostrings', update = FALSE, ask = FALSE); \
    if (!requireNamespace('Biostrings', quietly = TRUE)) quit(status = 1)"

RUN MAKEFLAGS="-j$(nproc)" R -e "\
    options(Ncpus = parallel::detectCores()); \
    install.packages('readxl', repos='https://cloud.r-project.org'); \
    if (!requireNamespace('readxl', quietly = TRUE)) quit(status = 1)"

# phylotools was archived on CRAN (2025-06-13, unresolved check issues) so
# install.packages() can no longer find it there. It's still installable
# straight from its CRAN source mirror on GitHub via remotes -- it only
# depends on `ape`, which is already pulled in as a dependency here.
RUN MAKEFLAGS="-j$(nproc)" R -e "\
    options(Ncpus = parallel::detectCores()); \
    install.packages('remotes', repos='https://cloud.r-project.org'); \
    remotes::install_github('cran/phylotools'); \
    if (!requireNamespace('phylotools', quietly = TRUE)) quit(status = 1)"

# -----------------------------------------------------------------------
# Sanity check at build time: fail the build early if anything is missing
# -----------------------------------------------------------------------
RUN R -e "library(tidyverse); library(Biostrings); library(phylotools); library(readxl)" \
    && mafft --version

# Nextflow runs containers non-interactively and expects a shell available;
# no ENTRYPOINT/CMD override needed, default bash is fine.
WORKDIR /data
