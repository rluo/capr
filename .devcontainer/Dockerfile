# =========================================
# Base: Ubuntu 24.04 (Noble)
# =========================================
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8

# 1) OS deps (+ build toolchain + libs R packages commonly need)
RUN apt-get update && apt-get install -y --no-install-recommends \
    locales ca-certificates wget git curl openssh-client \
    build-essential gfortran pkg-config \
    # R itself (from Ubuntu repo; version-agnostic as requested)
    r-base r-base-dev r-recommended \
    # C/C++/Fortran-compiled CRAN pkg deps
    libxml2-dev libcurl4-openssl-dev libssl-dev libgit2-dev libicu-dev \
    libopenblas-dev liblapack-dev libarmadillo-dev libfontconfig1-dev \
    # helpful for graphics/text shaping if pulled indirectly
    libpng-dev libjpeg-dev libtiff5-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev \
    && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen \
    && locale-gen \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 2) Install CRAN packages (non-interactive, parallel)
#    Feel free to add/remove pkgs in the vector below.
RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=parallel::detectCores()); \
    install.packages(c( \
    "Rcpp","RcppArmadillo","xml2","curl","openssl","git2r", \
    "languageserver","testthat","roxygen2","rversions","urlchecker", \
    "styler","lintr","devtools" \
    ), dependencies=TRUE)'

RUN git config --global user.email "xi.rossi.luo@gmail.com" && \
    git config --global user.name "Xi (Rossi) LUO"

# 3) Default workdir (mount your project here in VS Code)
WORKDIR /capr

# 4) Default command
# CMD ["R"]
