#########################################
### DOCKER SINGLE-CELL IMAGE PLATFORM ###
### MULTI-STAGE BUILD (SIZE OPTIMIZED) ###
#########################################
# Restructured path version: build context at repo root
# - R installers under docker/base/R/
# - Python requirements under docker/requirements/
# - Helper scripts under scripts/

#####################################
### STAGE 1: BUILDER (compile everything)
#####################################
FROM ubuntu:22.04 AS builder

LABEL maintainers="Anton Zhelonkin (anton.bioinf.md@gmail.com)"

####################################################
### Arguments and environment variables
####################################################
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# R version
ARG R_VERSION_MAJOR=4
ARG R_VERSION_MINOR=5
ARG R_VERSION_PATCH=3
ARG CRAN_MIRROR=https://cran.rstudio.com/

ENV R_VERSION_MAJOR=$R_VERSION_MAJOR \
    R_VERSION_MINOR=$R_VERSION_MINOR \
    R_VERSION_PATCH=$R_VERSION_PATCH \
    CONFIGURE_OPTIONS="--with-cairo --with-jpeglib --enable-R-shlib --with-blas --with-lapack" \
    CRAN_MIRROR=$CRAN_MIRROR

# GitHub token (build arg only)
ARG GITHUB_PAT

# User identity (for permission matching)
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=devuser
ARG GROUP=devgroup

####################################################
### System dependencies (builder stage)
####################################################
RUN apt-get update && apt-get install -y --no-install-recommends \
    tzdata \
    locales \
    build-essential g++ gfortran \
    git wget curl ca-certificates gnupg software-properties-common \
    # R build dependencies
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev \
    libuv1-dev libwebp-dev libnlopt-dev \
    libreadline-dev libbz2-dev liblzma-dev libicu-dev libpcre3-dev zlib1g-dev \
    libblas-dev liblapack-dev gfortran \
    # Additional scientific libraries (for R packages)
    libhdf5-dev \
    libgsl-dev \
    # gmp/glpk: gmp -> ClusterR -> mbkmeans, and igraph/leidenAlg -> rliger.
    # Their absence failed those packages silently.
    libgmp-dev libglpk-dev \
    # ImageMagick (for magick R package)
    libmagick++-dev \
    # X11 dependencies (required for R)
    libx11-dev libxt-dev x11proto-core-dev \
    # CLI tool build dependencies
    libbz2-dev liblzma-dev libncurses5-dev libcurl4-openssl-dev autoconf automake \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ >/etc/timezone \
    && locale-gen en_US.UTF-8 && update-locale LANG=en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

# Python 3.11 via deadsnakes PPA (Ubuntu 22.04 ships only 3.10; modern PyData
# packages — lightning, scvi-tools 1.2+, jax, flax — require >=3.11).
RUN add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-dev python3.11-venv python3.11-distutils \
    && rm -rf /var/lib/apt/lists/* && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11

####################################################
### Build R from source
####################################################
WORKDIR /tmp/r-build
RUN wget https://cran.rstudio.com/src/base/R-${R_VERSION_MAJOR}/R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz && \
    tar zxf R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz && \
    cd R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH} && \
    ./configure ${CONFIGURE_OPTIONS} && \
    make -j$(nproc) && \
    make install && \
    cd / && rm -rf /tmp/r-build

# Site-wide CRAN mirror
RUN mkdir -p /etc/R && \
    printf 'options(repos=c(CRAN="https://cloud.r-project.org"))\n' > /etc/R/Rprofile.site

# cmake >= 3.24, required by RcppPlanc (-> rliger). Ubuntu 22.04 apt ships
# 3.22.1, which RcppPlanc's CMakeLists rejects outright, so take it from pip.
RUN python3.11 -m pip install --no-cache-dir cmake && \
    ln -sf "$(python3.11 -c 'import cmake,os;print(os.path.join(os.path.dirname(cmake.__file__),"data","bin","cmake"))')" \
           /usr/local/bin/cmake && \
    cmake --version

# Packages declaring `CXX_STD = CXX11` (motifmatchr, chromVAR, Rfast, mbkmeans,
# rliger) compile at gnu++11, which current RcppArmadillo rejects outright.
# Makevars.site is read after Makeconf, so this raises the floor to gnu++17.
# Must stay in the builder stage, before any package installs.
RUN printf '%s\n' \
      'CXX11STD = -std=gnu++17' \
      'CXX14STD = -std=gnu++17' \
      > /usr/local/lib/R/etc/Makevars.site

####################################################
### Install R packages (via renv or initial core script)
####################################################
WORKDIR /opt/settings
COPY docker/base/R/install_core.R /opt/settings/install_R_packages.R
COPY docker/base/R/install_renv_project.R /opt/settings/
COPY docker/base/R/install_httpgd.R /opt/settings/

# Install R packages with aggressive cleanup
# Use --mount=type=secret to keep GITHUB_PAT out of layers
RUN --mount=type=secret,id=github_pat \
    if [ -f /run/secrets/github_pat ]; then \
        export GITHUB_PAT=$(cat /run/secrets/github_pat); \
    fi && \
    R -q --vanilla -f /opt/settings/install_renv_project.R && \
    R -q --vanilla -f /opt/settings/install_httpgd.R && \
    rm -rf ~/.cache/R/renv ~/.local/share/renv/cache /root/.cache/R /tmp/Rtmp* /tmp/downloaded_packages /tmp/Rd2pdf* /tmp/file*

# Generate manifest
RUN R -q -e 'pk <- as.data.frame(installed.packages()[, c("Package","Version","LibPath")]); \
             write.csv(pk, "/opt/settings/R-packages-manifest.csv", row.names = FALSE)'

####################################################
### Install TinyTeX (lightweight)
####################################################
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh && \
    /root/.TinyTeX/bin/*/tlmgr install collection-latex collection-fontsrecommended \
    collection-latexextra fancyvrb framed titling ae parskip xcolor && \
    /root/.TinyTeX/bin/*/tlmgr path add && \
    rm -rf /tmp/*

####################################################
### Build Python virtual environment
####################################################
COPY docker/requirements /opt/environments

ENV PIP_NO_CACHE_DIR=1
RUN python3.11 -m venv /opt/venvs/base && \
    /opt/venvs/base/bin/python -m pip install --upgrade pip setuptools wheel && \
    /opt/venvs/base/bin/python -m pip install --no-cache-dir -r /opt/environments/base.txt && \
    /opt/venvs/base/bin/python /opt/environments/verify_base.py && \
    rm -rf ~/.cache/pip /root/.cache/pip /tmp/*

####################################################
### Build bioinformatics CLI tools
####################################################
WORKDIR /build

# HTSlib, samtools, bcftools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xf htslib-1.21.tar.bz2 && cd htslib-1.21 && ./configure && make -j$(nproc) && make install && cd .. && \
    tar -xf bcftools-1.21.tar.bz2 && cd bcftools-1.21 && make -j$(nproc) && make install && cd .. && \
    tar -xf samtools-1.21.tar.bz2 && cd samtools-1.21 && make -j$(nproc) && make install && cd .. && \
    rm -rf /build/*

# BEDtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xf bedtools-2.31.1.tar.gz && cd bedtools2 && make -j$(nproc) && make install && \
    cd / && rm -rf /build/*

# scIBD (scATAC doublet detection, installed into the base venv).
# Upstream ships a FULL frozen conda env as its requirements.txt — hard `==`
# pins on numpy/pandas/scanpy/anndata/scikit-bio/... Installing it wholesale
# DOWNGRADES and breaks the carefully pinned base stack (scanpy's import chain
# fails: anndata -> dask -> xarray vs the downgraded pandas). So install scIBD
# itself with --no-deps and add only its genuine extra deps not already in
# base.txt, constrained by base.txt so the core stack is never downgraded.
# scIBD runs fine against the newer base scanpy/anndata.
RUN git clone https://github.com/Ying-Lab/scIBD && \
    /opt/venvs/base/bin/pip install --no-cache-dir --no-deps ./scIBD && \
    /opt/venvs/base/bin/pip install --no-cache-dir -c /opt/environments/base.txt \
        scikit-bio biom-format get-annotations annoy && \
    cd / && rm -rf /build/scIBD /root/.cache/pip /tmp/*

####################################################
### Final cleanup in builder stage
####################################################
RUN rm -rf /build /tmp/* /var/tmp/* /var/cache/apt/* /var/lib/apt/lists/*


#####################################
### STAGE 2: RUNTIME (clean base, copy artifacts)
#####################################
FROM ubuntu:22.04 AS runtime

# TODO(cache, deferred 2026-07): the runtime stage runs 4-5 separate
#   `apt-get update && install` layers (base deps, deadsnakes, tmux build-deps,
#   quarto, node). Each refreshes the apt index and adds a layer. Consolidating
#   the apt installs (single update, grouped installs where install order
#   allows) would trim layers and build time. Deferred: low risk but touches
#   the whole runtime stage, so batch it with the next structural change.

# IMAGE_VERSION is fed from the VERSION file by scripts/build.sh (--build-arg).
# Kept version-agnostic here so the label never drifts on a bump.
ARG IMAGE_VERSION=dev
LABEL maintainers="Anton Zhelonkin (anton.bioinf.md@gmail.com)"
LABEL version="${IMAGE_VERSION}"
LABEL description="Size-optimized single-cell bioinformatics environment (multi-stage build)"

####################################################
### Arguments and environment variables (runtime)
####################################################
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

ARG R_VERSION_MAJOR=4
ARG R_VERSION_MINOR=5
ARG R_VERSION_PATCH=3
ARG VENV_NAME=python3.11

ENV R_VERSION_MAJOR=$R_VERSION_MAJOR \
    R_VERSION_MINOR=$R_VERSION_MINOR \
    R_VERSION_PATCH=$R_VERSION_PATCH \
    VENV_NAME=$VENV_NAME \
    CRAN_MIRROR=https://cloud.r-project.org \
    IMAGE_FLAVOR=dev-core \
    PATH="/opt/venvs/base/bin:/root/.TinyTeX/bin/x86_64-linux:/usr/local/bin:${PATH}" \
    PIP_NO_CACHE_DIR=1 \
    VIRTUALENV_WRAPPER=/opt/venvs

####################################################
### User definition (runtime)
####################################################
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=devuser
ARG GROUP=devgroup

RUN apt-get update && apt-get install -y --no-install-recommends \
    tzdata \
    locales sudo \
    # Keep build-essential and dev libraries for RUNTIME package installation
    build-essential g++ gfortran \
    git wget curl ca-certificates \
    # R runtime dependencies
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev \
    libuv1-dev libwebp-dev libnlopt-dev \
    libreadline-dev libbz2-dev liblzma-dev libicu-dev libpcre3-dev zlib1g-dev \
    libblas-dev liblapack-dev \
    # Additional scientific libraries (for runtime R package compilation)
    libhdf5-dev \
    libgsl-dev \
    # Must mirror the builder: igraph links libglpk and gmp links libgmp, so
    # without these the runtime image has the .so files' owners missing and
    # igraph fails to load, taking SeuratData/GeneNMF/leidenAlg with it.
    libgmp-dev libglpk-dev \
    # ImageMagick (for magick R package)
    libmagick++-dev \
    # X11 dependencies (required for R)
    libx11-dev libxt-dev x11proto-core-dev \
    # Tooling for adding deadsnakes PPA below
    gnupg software-properties-common \
    # CLI tool runtime dependencies
    libncurses5-dev libcurl4-openssl-dev \
    # Utilities (ripgrep -> rg; fd-find -> fdfind, symlinked to fd below)
    # jq: JSON CLI required by SciAgent-toolkit's claude_settings.sh to backfill
    #     statusLine + gold defaults into .claude/settings.json (soft-fails without it).
    # pandoc: needed by htmlwidgets::saveWidget(selfcontained=TRUE) and Rmd/knitr
    #     rendering; previously best-effort apt-installed at runtime by install_deps.R.
    # imagemagick (convert/identify): image.nvim's magick_cli processor for
    #     inline figures. chafa: sixel/ansi fallback. tldr: quick command help.
    #     See docs/plans/inline-images-nvim.md (dev-env repo).
    # NOTE: tmux is intentionally NOT apt-installed here — Ubuntu 22.04 ships
    #     3.2a, which lacks `allow-passthrough` (the kitty-graphics switch). A
    #     pinned modern tmux is built from source below into /usr/local/bin.
    vim nano less htop ripgrep fd-find jq pandoc \
    imagemagick chafa tldr \
    # bubblewrap (bwrap): the sandbox backend for the Codex CLI's
    #   workspace-write mode. Without it, sandboxed exec fails outright. Note it
    #   still needs unprivileged user namespaces at RUN time — pair with
    #   `--security-opt seccomp=unconfined` (see docs/ai-integration.md).
    bubblewrap && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ >/etc/timezone && \
    locale-gen en_US.UTF-8 && update-locale LANG=en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

# Python 3.11 runtime (matches builder venv interpreter)
RUN add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-dev python3.11-venv python3.11-distutils \
    && rm -rf /var/lib/apt/lists/*

# Expose Debian's fd-find binary (installed as `fdfind`) under its canonical
# `fd` name, system-wide on PATH (/usr/local/bin is already on PATH above).
RUN ln -sf "$(command -v fdfind)" /usr/local/bin/fd

####################################################
### tmux >= 3.3 (source build) — kitty inline-image passthrough
####################################################
# Ubuntu 22.04 apt only ships tmux 3.2a, which LACKS `allow-passthrough` — the
# switch that lets the kitty graphics protocol (image.nvim) tunnel through tmux
# to the outer terminal. Build a pinned modern tmux from source into
# /usr/local/bin (first on PATH) so inline figures work inside the container.
# See docs/plans/inline-images-nvim.md (dev-env repo). ncurses headers come from
# libncurses5-dev (installed above); libevent + bison are pulled just for this.
# Pinned to the latest official release (matches the workstation host tmux), so
# the container and host multiplexers are identical. Any >= 3.3 has passthrough.
#
# TODO(size, deferred 2026-07): this RUN leaves `bison` and `libevent-dev` in
#   the runtime image (build-only). Purging them shaves a little bloat, but
#   `libevent-dev` pulls the `libevent-2.1-7` runtime lib that tmux links
#   against — naive `apt-get purge` + autoremove would break tmux. Correct fix:
#   purge `bison` + `libevent-dev` but explicitly reinstall `libevent-2.1-7`
#   in the same layer. Skipped this pass to avoid a tmux-breakage risk.
ARG TMUX_VERSION=3.7b
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends libevent-dev bison; \
    cd /tmp; \
    curl -fsSL -o tmux.tar.gz "https://github.com/tmux/tmux/releases/download/${TMUX_VERSION}/tmux-${TMUX_VERSION}.tar.gz"; \
    tar xzf tmux.tar.gz; \
    cd "tmux-${TMUX_VERSION}"; \
    ./configure --prefix=/usr/local; \
    make -j"$(nproc)"; \
    make install; \
    cd /; rm -rf /tmp/tmux*; \
    rm -rf /var/lib/apt/lists/*; \
    tmux -V

####################################################
### Quarto CLI (native .qmd render for scientific notebooks)
####################################################
# Quarto ships as a standalone .deb from GitHub (bundles its own pandoc + engine);
# it is NOT in the Ubuntu archives, so `apt-get install quarto` never resolves.
# Pinned for reproducibility. Enables `quarto render`; the apt `pandoc` above still
# provides a plain `pandoc` on PATH for rmarkdown/knitr + htmlwidgets::saveWidget.
ARG QUARTO_VERSION=1.9.38
RUN set -eux; \
    arch="$(dpkg --print-architecture)"; \
    case "$arch" in amd64) qarch=amd64 ;; arm64) qarch=arm64 ;; *) echo "unsupported arch: $arch" >&2; exit 1 ;; esac; \
    curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${qarch}.deb"; \
    apt-get update; \
    apt-get install -y --no-install-recommends /tmp/quarto.deb; \
    rm -f /tmp/quarto.deb; \
    rm -rf /var/lib/apt/lists/*; \
    quarto --version

####################################################
### MALLET (LDA topic modelling backend for pycisTopic)
####################################################
# LDA backend for pycisTopic's run_cgs_models_mallet(). Java-only, and the
# image's sole Java consumer, hence a headless JRE. The stock 1g heap OOMs on
# realistic models, so MALLET_MEMORY is raised; mallet must stay on PATH.
ARG MALLET_VERSION=2.0.8
ENV MALLET_HOME=/opt/mallet \
    MALLET_MEMORY=16g
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends openjdk-17-jre-headless; \
    rm -rf /var/lib/apt/lists/*; \
    cd /tmp; \
    curl -fsSL -o mallet.tar.gz \
      "https://mallet.cs.umass.edu/dist/mallet-${MALLET_VERSION}.tar.gz"; \
    mkdir -p "${MALLET_HOME}"; \
    tar xzf mallet.tar.gz -C "${MALLET_HOME}" --strip-components=1; \
    rm -f mallet.tar.gz; \
    chmod +x "${MALLET_HOME}/bin/mallet"; \
    ln -sf "${MALLET_HOME}/bin/mallet" /usr/local/bin/mallet; \
    java -version; \
    mallet --help >/dev/null 2>&1 || true; \
    test -f "${MALLET_HOME}/class/cc/mallet/topics/ParallelTopicModel.class" \
      || test -f "${MALLET_HOME}/lib/mallet.jar"

####################################################
### AI prerequisites (Node.js + uv)
####################################################
# These are infrastructure prerequisites only. The image carries NO AI tooling
# itself (no Claude/Gemini CLI, no MCP servers, no ToolUniverse/Serena/PAL).
# All AI tooling is installed at runtime per-project by:
#   01_modules/SciAgent-toolkit/scripts/setup-ai.sh
# which depends on `node`/`npm`/`npx` and `uv`/`uvx` being available.

# Node.js 20 LTS from NodeSource (provides node, npm, npx)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get update && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/* && \
    node --version && npm --version && npx --version

# uv / uvx (Astral) — fast Python package/runtime manager
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    mv /root/.local/bin/uv /usr/local/bin/ && \
    mv /root/.local/bin/uvx /usr/local/bin/ && \
    uv --version

# Create user with proper permissions
RUN set -eux; \
    SAN_GROUP="$(printf '%s' "${GROUP}" | tr '[:space:]' '_')"; \
    if getent group "${GROUP_ID}" >/dev/null; then \
        GROUP_NAME="$(getent group "${GROUP_ID}" | cut -d: -f1)"; \
    else \
        groupadd -g "${GROUP_ID}" "${SAN_GROUP}" || groupadd -g "${GROUP_ID}" "grp_${GROUP_ID}"; \
        GROUP_NAME="$(getent group "${GROUP_ID}" | cut -d: -f1)"; \
    fi; \
    if id -u "${USER_ID}" >/dev/null 2>&1; then \
        USER_NAME="$(getent passwd "${USER_ID}" | cut -d: -f1)"; \
    else \
        useradd -m -s /bin/bash -u "${USER_ID}" -g "${GROUP_ID}" "${USER}"; \
        USER_NAME="${USER}"; \
    fi; \
    mkdir -p /home/${USER_NAME}/reference_genome /home/${USER_NAME}/data; \
    chown -R ${USER_ID}:${GROUP_ID} /home/${USER_NAME}; \
    echo "alias ll='ls -la --color=auto'" >> /home/${USER_NAME}/.bashrc; \
    mkdir -p /etc/sudoers.d; \
    echo "${USER_NAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-${USER_NAME}; \
    chmod 440 /etc/sudoers.d/90-${USER_NAME}

####################################################
### Copy compiled artifacts from builder
####################################################

# R installation
COPY --from=builder /usr/local/lib/R /usr/local/lib/R
COPY --from=builder /usr/local/bin/R /usr/local/bin/R
COPY --from=builder /usr/local/bin/Rscript /usr/local/bin/Rscript
COPY --from=builder /etc/R /etc/R

# TinyTeX
COPY --from=builder /root/.TinyTeX /root/.TinyTeX

# Python venv (base)
COPY --from=builder /opt/venvs/base /opt/venvs/base

####################################################
### Interactive Python widgets (jupyter-scatter + live kernel)
####################################################
# jupyter-scatter (a.k.a. jscatter) renders interactive scatterplots in the
# notebook via anywidget/ipywidgets; those require a LIVE Jupyter kernel, hence
# ipykernel. We install `jupyter-scatter` — the lean meta-package — NOT the raw
# `jscatter` distribution (which drags in MDAnalysis/nglview/pdb2pqr). Versions
# pinned to what the running STING-JR base venv already validated. This runs
# AFTER the venv is COPY'd from the builder, so /opt/venvs/base/bin/pip is valid.
#
# The kernelspec is registered with --prefix (into the venv's share/jupyter),
# NEVER --user: the runtime user is `devuser` (uid 1000) and cannot see kernels
# installed under /root/.local. --prefix keeps the spec inside the venv, visible
# to every user (fixes VS Code's "Running cells requires ipykernel" prompt).
RUN /opt/venvs/base/bin/pip install --no-cache-dir \
        "jupyter-scatter==1.0.1" \
        "ipykernel==7.3.0" && \
    /opt/venvs/base/bin/python -m ipykernel install \
        --prefix /opt/venvs/base \
        --name python311-scagent \
        --display-name "Python 3.11 (SciAgent)" && \
    rm -rf /root/.cache/pip ~/.cache/pip /tmp/* 2>/dev/null || true

####################################################
### Register the R Jupyter kernel (IRkernel)
####################################################
# install_core.R installs the IRkernel *package* but never registers a
# kernelspec, so jupyter_client (and molten.nvim) cannot see an R kernel.
# Register it with prefix= into the venv's share/jupyter — the same location
# and reasoning as the Python kernel above: a --user spec lands in /root/.local
# and is invisible to devuser (uid 1000). See docs/plans/inline-images-nvim.md.
RUN R -q --vanilla -e 'IRkernel::installspec(prefix = "/opt/venvs/base", name = "ir", displayname = "R (SciAgent)")' && \
    rm -rf /tmp/* 2>/dev/null || true

# Bioinformatics tools (samtools, bcftools, bedtools)
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=builder /usr/local/bin/bedtools /usr/local/bin/bedtools
COPY --from=builder /usr/local/lib/libhts* /usr/local/lib/
COPY --from=builder /usr/local/include/htslib /usr/local/include/htslib

# R settings and scripts
COPY --from=builder /opt/settings /opt/settings

# Environment files for layered venvs
COPY docker/requirements /opt/environments

# Layered venv helper script
COPY scripts/create_layered_venv.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/create_layered_venv.sh

# NOTE: SciAgent-toolkit is intentionally NOT copied into the image.
# It is added per-project as a git submodule under 01_modules/SciAgent-toolkit/
# and its setup-ai.sh runs inside the container at project-init time.

# Update shared library cache
RUN ldconfig

####################################################
### User configuration (R profiles, venv wrappers)
####################################################

# Copy .Rprofile for both root and user (source file currently under .devcontainer)
COPY .devcontainer/.Rprofile /root/
RUN set -eux; \
    user_home="$(getent passwd ${USER_ID} | cut -d: -f6)"; \
    cp /root/.Rprofile "${user_home}/.Rprofile"; \
    chown ${USER_ID}:${GROUP_ID} "${user_home}/.Rprofile"; \
    mkdir -p "${user_home}/.vscode-R" /root/.vscode-R; \
    printf '# VS Code R init (placeholder)\n' > "${user_home}/.vscode-R/init.R"; \
    printf '# VS Code R init (placeholder)\n' > /root/.vscode-R/init.R; \
    chown -R ${USER_ID}:${GROUP_ID} "${user_home}/.vscode-R"; \
    # Create user R library for runtime installs
    mkdir -p "${user_home}/R/x86_64-pc-linux-gnu-library/${R_VERSION_MAJOR}.${R_VERSION_MINOR}"; \
    chown -R ${USER_ID}:${GROUP_ID} "${user_home}/R"

####################################################
### Python venv wrappers
####################################################
RUN printf '%s\n' \
  '#!/usr/bin/env bash' \
  'source /opt/venvs/base/bin/activate && exec "$@"' \
  > /usr/local/bin/py-base && chmod +x /usr/local/bin/py-base && \
  printf '%s\n' \
  '# usage: usepy base|squid|atac|comms|scenic' \
  '# Switches to layered venv (create if needed)' \
  'usepy(){ ' \
  '  case "$1" in ' \
  '    base) . /opt/venvs/base/bin/activate ;;' \
  '    squid|atac|comms) ' \
  '      if [ ! -d "/opt/venvs/$1" ]; then ' \
  '        echo "Creating layered venv: $1"; ' \
  '        python3.11 -m venv --system-site-packages /opt/venvs/$1; ' \
  '        if [ -f "/opt/environments/${1}.txt" ]; then ' \
  '          /opt/venvs/$1/bin/pip install --no-cache-dir -r /opt/environments/${1}.txt; ' \
  '        fi; ' \
  '      fi; ' \
  '      . /opt/venvs/$1/bin/activate ;;' \
  '    scenic) ' \
  '      # ISOLATED, not layered: pycisTopic pins pandas==1.5 and must not' \
  '      # shadow base pandas 2.2.3 for scanpy/anndata/muon. See scenic.txt.' \
  '      if [ ! -d "/opt/venvs/scenic" ]; then ' \
  '        echo "Creating ISOLATED venv: scenic (no base inheritance; this takes a while)"; ' \
  '        python3.11 -m venv /opt/venvs/scenic; ' \
  '        /opt/venvs/scenic/bin/pip install --no-cache-dir --upgrade pip setuptools wheel; ' \
  '        if [ -f /opt/environments/scenic.txt ]; then ' \
  '          /opt/venvs/scenic/bin/pip install --no-cache-dir -r /opt/environments/scenic.txt; ' \
  '        fi; ' \
  '      fi; ' \
  '      . /opt/venvs/scenic/bin/activate ;;' \
  '    *) echo "usepy base|squid|atac|comms|scenic"; return 1;;' \
  '  esac; ' \
  '}' \
  '# Put user-local binaries (e.g. dev-env nvim/fzf in ~/.local/bin) on PATH' \
  'case ":$PATH:" in *":$HOME/.local/bin:"*) ;; *) export PATH="$HOME/.local/bin:$PATH";; esac' \
  >> /etc/bash.bashrc

####################################################
### R session wrappers
####################################################
RUN printf '%s\n' \
  '#!/usr/bin/env bash' \
  'unset USE_ARCHR' \
  'exec /opt/venvs/base/bin/radian "$@"' \
  > /usr/local/bin/r-base && chmod +x /usr/local/bin/r-base && \
  printf '%s\n' \
  '#!/usr/bin/env bash' \
  'export USE_ARCHR=1' \
  'exec /opt/venvs/base/bin/radian "$@"' \
  > /usr/local/bin/r-archr && chmod +x /usr/local/bin/r-archr

####################################################
### Final permissions and workspace
####################################################
RUN mkdir -p /workspaces/project && \
    chown -R ${USER_ID}:${GROUP_ID} /workspaces/project

WORKDIR /workspaces/project
USER ${USER_ID}:${GROUP_ID}

CMD ["/bin/bash"]

