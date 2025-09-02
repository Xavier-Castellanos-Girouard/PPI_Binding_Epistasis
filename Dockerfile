# Use Ubuntu 20.04 as the base image
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

# Install Python 3.9 and required system dependencies
RUN apt-get update && apt-get install -y \
    python3.9 \
    python3.9-dev \
    python3.9-distutils \
    python3.9-venv \
    ca-certificates \
    wget \
    gnupg \
    libxml2-dev \
    libxslt1-dev \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.9 as the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 1

# Ensure Python 3.9 has pip installed
RUN python3.9 -m ensurepip --upgrade || curl -sS https://bootstrap.pypa.io/get-pip.py | python3.9

# Upgrade pip
RUN python3.9 -m pip install --upgrade pip

# Install required Python packages
RUN python3.9 -m pip install --no-cache-dir \
    numpy==1.24.1 \
    matplotlib==3.5.2 \
    networkx==2.8.5 \
    seaborn==0.11.2 \
    pandas==1.4.3 \
    tqdm==4.64.0 \
    scipy==1.8.1 \
    requests==2.28.1 \
    lxml==4.5.0 \
    beautifulsoup4==4.8.2 \
    openpyxl==3.0.10

# Install necessary packages including ca-certificates, wget, Java, and build tools for Cluster 3
RUN apt-get update && apt-get install -y \
    openjdk-11-jre \
    libx11-6 \
    libxt6 \
    libxext6 \
    build-essential \
    automake \
    autoconf \
    m4 \
    perl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxslt1-dev \
    libgmp-dev \
    libmpfr-dev \
    cmake \
    clustalw \
    && rm -rf /var/lib/apt/lists/*


# Add the CRAN repository and GPG key for R 4.2
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /usr/share/keyrings/cran-archive-keyring.gpg > /dev/null \
    && echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" | tee /etc/apt/sources.list.d/cran.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       r-base=4.2.3-1.2004.0 \
       r-base-dev=4.2.3-1.2004.0 \
       r-recommended=4.2.3-1.2004.0 \
    && rm -rf /var/lib/apt/lists/*

# Ensure R is installed correctly
RUN R --version

# Install essential R packages
RUN R -e "install.packages(c('dplyr', 'tidyr', 'ggplot2', 'stringr', 'seqinr', 'data.table', 'readxl', 'ggpubr', 'cowplot', 'ggsignif', 'viridis', 'ggstatsplot'), repos='https://cloud.r-project.org/')"

# Install BiocManager and org.Sc.sgd.db from Bioconductor
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('curl', 'openssl', 'httr', 'AnnotationDbi', 'org.Sc.sgd.db'), update=FALSE, ask=FALSE)"
