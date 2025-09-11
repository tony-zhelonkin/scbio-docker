#########################################
### DOCKER SINGLE-CELL IMAGE PLATFORM ###
#########################################
# Step 1 - Import base image
#############################
FROM ubuntu:22.04 
# Label
LABEL maintainers="Anton Zhelonkin (anton.bioinf.md@gmail.com)"

####################################################
### Step 2: Set arguments and environment variables
####################################################

# Set noninteractive mode
ARG DEBIAN_FRONTEND=noninteractive

# Timezone definition
ENV TZ=America/Chicago LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ >/etc/timezone && \
    apt-get update && apt-get install -y --no-install-recommends locales && \
    locale-gen en_US.UTF-8 && update-locale LANG=en_US.UTF-8

# Tools definition
# Python environment
ARG VENV_NAME=python3.10
# Define environment variables
ENV VENV_NAME=$VENV_NAME
# Keep all Python environment definitions inside the image
COPY .environments /opt/environments


# R version pinned by default
ARG R_VERSION_MAJOR=4
ARG R_VERSION_MINOR=4
ARG R_VERSION_PATCH=2
ARG CRAN_MIRROR=https://cran.rstudio.com/
ARG QUARTO_VER="1.6.39"

# Define R environment variables
ENV R_VERSION_MAJOR=$R_VERSION_MAJOR
ENV R_VERSION_MINOR=$R_VERSION_MINOR
ENV R_VERSION_PATCH=$R_VERSION_PATCH
ENV QUARTO_VER=$QUARTO_VER
ENV CONFIGURE_OPTIONS="--with-cairo --with-jpeglib --enable-R-shlib --with-blas --with-lapack"
ENV CRAN_MIRROR=$CRAN_MIRROR

# GitHub API key
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT


############################
### Step 3: User definition
############################
# Pass host identity at build if you like:
#   --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) \
#   --build-arg USER=$USER --build-arg GROUP=$(id -gn)
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=devuser
ARG GROUP=devgroup

# Create group/user matching host's UID/GID (robust to spaces in GROUP)
RUN set -eux; \
    SAN_GROUP="$(printf '%s' "${GROUP}" | tr '[:space:]' '_')"; \
    if getent group "${GROUP_ID}" >/dev/null; then \
        GROUP_NAME="$(getent group "${GROUP_ID}" | cut -d: -f1)"; \
    else \
        # If creation with provided name fails (e.g., illegal chars), fall back to grp_<gid>
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
    # Use numeric IDs here to avoid name parsing issues
    chown -R ${USER_ID}:${GROUP_ID} /home/${USER_NAME}; \
    echo "alias ll='ls -la --color=auto'" >> /home/${USER_NAME}/.bashrc; \
    # ensure sudoers.d exists, then add NOPASSWD rule for this user
    mkdir -p /etc/sudoers.d; \
    echo "${USER_NAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-${USER_NAME}; \
    chmod 440 /etc/sudoers.d/90-${USER_NAME}

###############################
### Step 4: System dependencies
###############################

RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build essentials and compilers
    build-essential \
    g++ \
    gfortran \
    make \
    cmake \
    automake \
    \
    # Python ecosystem
    python3.10 \
    python3.10-dev \
    python3.10-venv \
    python3-pip \
    python3-wheel \
    python3-setuptools \
    python3-launchpadlib \
    \
    # Development libraries
    perl \
    autoconf \
    libtool \
    pkg-config \
    libboost-all-dev \
    libblosc-dev \
    liblzo2-dev \
    libhdf5-serial-dev \
    libreadline-dev \
    readline-common \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libgit2-dev \
    libglpk-dev \
    libgsl-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    liblapack-dev \
    liblzma-dev \
    libmagick++-dev \
    libncurses-dev \
    libomp-dev \
    libopenblas-dev \
    libpcre2-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libtbb-dev \
    libzmq3-dev \
    libudunits2-dev \
    libgdal-dev \
    zlib1g-dev \
    libghc-zlib-dev \
    libbz2-dev \
    libzstd-dev \
    \
    # System utilities
    apt-utils \
    ca-certificates \
    curl \
    file \
    pandoc \
    less \
    bison \
    git \
    zip \
    gzip \
    gnupg \
    locales \
    openssh-server \
    openssl \
    binutils \
    parallel \
    screen \
    shellinabox \
    software-properties-common \
    sudo \
    supervisor \
    tmux \
    tzdata \
    unzip \
    vim \
    wget \
    ftp \
    xvfb \
    \
    # TeX support
    texinfo \
    texlive \
    texlive-fonts-extra \
    \
    # Java
    default-jdk \
    maven \
    openjdk-17-jre-headless \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


###########################
### Step 5: R installation
###########################
RUN wget https://cran.rstudio.com/src/base/R-${R_VERSION_MAJOR}/R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz && \
    tar zxvf R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz && \
    rm R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz

WORKDIR /R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}
RUN ./configure ${CONFIGURE_OPTIONS} && \
    make && \
    make install

RUN locale-gen en_US.UTF-8
WORKDIR /root
# Put settings & scripts somewhere stable
RUN mkdir -p /opt/settings
COPY .devcontainer/install_R_packages.R /opt/settings/install_R_packages.R
COPY .devcontainer/install_quarto.sh /opt/settings/install_quarto.sh
COPY .devcontainer/.Rprofile /root/
RUN mkdir -p /root/.vscode-R && echo '# VS Code R init (placeholder)' > /root/.vscode-R/init.R
RUN set -eux; \
    user_home="$(getent passwd ${USER_ID} | cut -d: -f6)"; \
    cp /root/.Rprofile "${user_home}/.Rprofile"; \
    chown ${USER_ID}:${GROUP_ID} "${user_home}/.Rprofile"; \
    mkdir -p "${user_home}/.vscode-R"; \
    printf '# VS Code R init (placeholder)\n' > "${user_home}/.vscode-R/init.R"; \
    chown -R ${USER_ID}:${GROUP_ID} "${user_home}/.vscode-R"

# If you have a lock pre-committed, copy it here; if not, this line is harmless.
# (You can add /opt/settings/renv.lock later and the restore step below will use it.)
# COPY renv.lock /opt/settings/renv.lock

# Install Quarto
RUN bash /opt/settings/install_quarto.sh $QUARTO_VER

# Install renv and then either:
#  - restore from lockfile if present, OR
#  - run your installer and snapshot to create /opt/settings/renv.lock
# RUN R -q -e 'install.packages("renv", repos="https://cloud.r-project.org")' && \
#     R -q -e 'if (file.exists("/opt/settings/renv.lock")) { \
#                renv::init(bare = TRUE); \
#                renv::restore(lockfile = "/opt/settings/renv.lock", prompt = FALSE); \
#              } else { \
#                source("/opt/settings/install_R_packages.R"); \
#                renv::init(bare = TRUE); \
#                renv::snapshot(lockfile = "/opt/settings/renv.lock", prompt = FALSE); \
#              }'
# was: R -q -e 'install.packages("renv", ...)'
# and:  R -q -e 'if (file.exists("/opt/settings/renv.lock")) { ... }'
RUN R -q --vanilla -e 'install.packages("renv", repos="https://cloud.r-project.org")' && \
    R -q --vanilla -e 'if (file.exists("/opt/settings/renv.lock")) { \
                         renv::init(bare = TRUE); \
                         renv::restore(lockfile = "/opt/settings/renv.lock", prompt = FALSE); \
                       } else { \
                         source("/opt/settings/install_R_packages.R"); \
                         renv::init(bare = TRUE); \
                         renv::snapshot(lockfile = "/opt/settings/renv.lock", prompt = FALSE); \
                       }'

# Human-readable manifest of installed packages (helps auditing outside renv)
RUN R -q -e 'pk <- as.data.frame(installed.packages()[, c("Package","Version","LibPath")]); \
             write.csv(pk, "/opt/settings/R-packages-manifest.csv", row.names = FALSE)'

#########################
### Step 6: Pip installs
#########################
#########################
### Step 6: Python envs
#########################
# Python is already installed via apt (3.10) above.
ENV PIP_NO_CACHE_DIR=1 \
    VIRTUALENV_WRAPPER=/opt/venvs

# Create both venvs
RUN python3 -m venv /opt/venvs/base && \
    python3 -m venv /opt/venvs/squid

# Upgrade pip & install requirements into each
RUN /opt/venvs/base/bin/python -m pip install --upgrade pip setuptools wheel && \
    /opt/venvs/base/bin/python -m pip install -r /opt/environments/base_requirements.txt

RUN /opt/venvs/squid/bin/python -m pip install --upgrade pip setuptools wheel && \
    /opt/venvs/squid/bin/python -m pip install -r /opt/environments/squid_requirements.txt

# Default to the "base" venv on PATH; easy to switch later
ENV PATH="/opt/venvs/base/bin:${PATH}"

# Small helper wrappers to switch envs quickly in interactive shells
RUN printf '%s\n' \
  '#!/usr/bin/env bash' \
  'source /opt/venvs/base/bin/activate && exec "$@"' \
  > /usr/local/bin/py-base && chmod +x /usr/local/bin/py-base && \
  printf '%s\n' \
  '#!/usr/bin/env bash' \
  'source /opt/venvs/squid/bin/activate && exec "$@"' \
  > /usr/local/bin/py-squid && chmod +x /usr/local/bin/py-squid && \
  printf '%s\n' \
  '# usage: usepy base|squid' \
  'usepy(){ case "$1" in base) . /opt/venvs/base/bin/activate ;; squid) . /opt/venvs/squid/bin/activate ;; *) echo "usepy base|squid"; return 1;; esac; }' \
  >> /etc/bash.bashrc

#   RUN pip3 install setuptools wheel 
#   RUN pip3 install numpy scipy pandas
#   RUN pip3 install anndata==0.10.9
#   RUN pip3 install "python-igraph==0.10.4"



#################################
### Step 7: Bioinformatics tools
#################################
# Creating the build directory
RUN mkdir -p /build && \
    chown -R ${USER_ID}:${GROUP_ID} /build && \
    chmod 755 /build
WORKDIR /build

### Bowtie2 (removed to slim image)

##############
### scIBD ### 
##############
RUN git clone https://github.com/Ying-Lab/scIBD
WORKDIR /build/scIBD
RUN pip install -r requirements.txt && \
    python3 setup.py install
WORKDIR /build


### fqtools (removed to slim image)

# Optionally drop build-time autotools if you won’t compile other C projects later:
# RUN apt-get purge -y --auto-remove automake autoconf libtool pkg-config && \
#     apt-get clean && rm -rf /var/lib/apt/lists/*

### Chromap (removed to slim image)


####################################
### HTSlib tools: Samtools, Bcftools
####################################
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xf htslib-1.21.tar.bz2

WORKDIR /build/htslib-1.21
RUN ./configure && make -j 4 && make install
WORKDIR /build
RUN tar -xf bcftools-1.21.tar.bz2
WORKDIR /build/bcftools-1.21
RUN make -j 4 && make install
WORKDIR /build
RUN tar -xf samtools-1.21.tar.bz2
WORKDIR /build/samtools-1.21
RUN make -j 4 && make install
WORKDIR /build

### STAR (removed to slim image)

### kallisto (removed to slim image)

### Salmon (removed to slim image)

#############
### BEDtools
#############
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xf bedtools-2.31.1.tar.gz && cd bedtools2 && \
    make -j 4 && make install

### Java tools (FastQC, Trimmomatic) removed to slim image

### Picard removed to slim image

### Trim Galore removed to slim image

### featureCounts removed to slim image

### BWA removed to slim image

### CHIP-seq tools (USeq, Sissrs) removed to slim image


###########################
### Cleaning up the build 
###########################
WORKDIR /build
RUN rm -rf /build/*

# Switch to non-root user by default
USER ${USER_ID}:${GROUP_ID}
ENV HOME=/home/${USER}
WORKDIR /home/${USER}