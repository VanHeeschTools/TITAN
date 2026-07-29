FROM rocker/shiny:4.4.2

# System libraries for compiled R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# remotes is the pinning mechanism; install it first
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Pinned to match local environment.
# upgrade='never' avoids touching already-installed packages; fresh base means
# all required dependencies are fetched from scratch without version conflicts.
RUN Rscript -e "\
  r  <- 'https://cloud.r-project.org'; \
  iv <- function(p, v) remotes::install_version(p, version=v, repos=r, upgrade='never'); \
  iv('shiny',        '1.14.0'); \
  iv('bslib',        '0.11.0'); \
  iv('DT',           '0.34.0'); \
  iv('plotly',       '4.12.0'); \
  iv('dplyr',        '1.2.1');  \
  iv('tidyr',        '1.3.2');  \
  iv('stringr',      '1.6.0');  \
  iv('ggplot2',      '4.0.3');  \
  iv('shinyWidgets', '0.9.1')"

# Unpinned utilities (no version constraint from local env)
RUN Rscript -e "install.packages(c('cowplot', 'data.table', 'jsonlite', 'yaml', 'processx', 'XML', 'BiocManager'), repos='https://cloud.r-project.org')"

# Bioconductor packages — rBLAST is a Bioconductor package, not CRAN
RUN Rscript -e "BiocManager::install(c('Biostrings', 'IRanges', 'rBLAST'), ask=FALSE, update=FALSE)"

COPY app/ /srv/titan/

ENV PATH="/srv/titan/bin:${PATH}"

RUN chown -R shiny:shiny /srv/titan

USER shiny

# Cloud Run injects $PORT (default 8080); local docker run falls back to 3838.
# Ref data and study data are provided via volume mounts at /srv/titan/ref and
# /srv/titan/data — either local bind mounts (dev) or GCS volumes (Cloud Run).
CMD ["Rscript", "-e", \
     "shiny::runApp('/srv/titan', host='0.0.0.0', port=as.integer(Sys.getenv('PORT', 3838)))"]
