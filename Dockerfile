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
RUN Rscript -e "install.packages(c('cowplot', 'data.table', 'jsonlite'), repos='https://cloud.r-project.org')"

# Serve TITAN at the root path
RUN printf 'run_as shiny;\n\nserver {\n  listen 3838;\n  location / {\n    app_dir /srv/titan;\n    log_dir /var/log/shiny-server;\n  }\n}\n' \
    > /etc/shiny-server/shiny-server.conf

COPY app/ /srv/titan/

RUN chown -R shiny:shiny /srv/titan

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=10s CMD \
    curl -sf http://localhost:3838/ || exit 1
