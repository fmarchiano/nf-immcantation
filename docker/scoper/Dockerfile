FROM r-base:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('remotes', repos='https://cran.r-project.org')"
RUN Rscript -e "remotes::install_version('alakazam', version='1.4.3', repos='https://cran.r-project.org')"
RUN Rscript -e "remotes::install_version('shazam',   version='1.3.2', repos='https://cran.r-project.org')"
RUN Rscript -e "remotes::install_version('scoper',   version='1.5.0', repos='https://cran.r-project.org')"
