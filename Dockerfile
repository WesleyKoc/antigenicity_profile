# Base image: Match your renv.lock R version (change 4.3.1 to match your target R version)
FROM rocker/rstudio:4.5.2

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    python3 \
    python3-pip \
    libncurses-dev \
    libdeflate-dev \
    libssl-dev \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    clustalo \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    libmagick++-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    muscle \
    && rm -rf /var/lib/apt/lists/*

# Set working directory for tools installation
WORKDIR /opt/tools

##############################
# Install samtools from source
##############################
RUN curl -L -O https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 \
    && tar -jxf samtools-1.22.tar.bz2 \
    && cd samtools-1.22 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd /opt/tools \
    && rm -rf samtools-1.22 samtools-1.22.tar.bz2

##############################
# Install bcftools from source
##############################
RUN curl -L -O https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -jxf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd /opt/tools \
    && rm -rf bcftools-1.22 bcftools-1.22.tar.bz2

##############################
# Install R Packages via renv
##############################
WORKDIR /srv/shiny-server/myapp

# Copy renv files first to leverage Docker layer caching
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/ renv/

# Bootstrap renv and restore exact package states from lockfile
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "renv::restore(prompt = FALSE)"

# Copy the rest of your application code (R/, data/, www/, result/, app.R)
COPY . .

# Expose RStudio and Shiny ports
EXPOSE 8787 3838

# Default command for rocker/rstudio images
CMD ["/init"]