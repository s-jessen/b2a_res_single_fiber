# Use rocker image with R 4.5.1
FROM rocker/r-ver:4.5.1

# Set environment variables to reduce prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install Renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Install system libraries (needed by many R packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    fonts-dejavu-core \
    libglpk-dev \
    libsuitesparse-dev \
    curl \
    git \
    wget \
    unzip \
    ca-certificates \
    build-essential && \
    # --- Install Source Sans 3 manually ---
    mkdir -p /usr/share/fonts/truetype/source-sans-pro && \
    cd /usr/share/fonts/truetype/source-sans-pro && \
    wget -q https://github.com/adobe-fonts/source-sans/releases/download/3.052R/OTF-source-sans-3.052R.zip && \
    unzip OTF-source-sans-3.052R.zip && \
    rm OTF-source-sans-3.052R.zip && \
    fc-cache -fv && \
    # --- Cleanup ---
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniconda
ENV PATH=/opt/conda/bin:$PATH
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda config --set always_yes yes && \
    /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --set channel_priority strict && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    /opt/conda/bin/conda update -n base conda

# Copy project into container ('data', 'data-raw' and 'R' folders ignored by .dockerignore
# because they are mounted at runtime)
WORKDIR /project
COPY . /project

# Restore the R environment
RUN R -e "renv::restore()"

# Install conda env using environment.yml
RUN conda env create -f environment.yml
