FROM r-base:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Debian's default CFLAGS include -Werror=format-security, which breaks the
# Bioconductor S4Vectors source build (transitive dep of alakazam). Override
# system-wide before any compiled package is installed.
RUN mkdir -p /root/.R && \
    printf 'CFLAGS=-g -O2 -Wno-format-security\nCXXFLAGS=-g -O2 -Wno-format-security\n' \
        > /root/.R/Makevars

RUN Rscript -e "install.packages(c('remotes','BiocManager'), repos='https://cran.r-project.org')"

# alakazam 1.4.3 imports Biostrings/GenomicAlignments/IRanges from Bioconductor;
# install those first via BiocManager so the CRAN install_version calls succeed.
RUN Rscript -e "BiocManager::install(c('Biostrings','GenomicAlignments','IRanges'), update=FALSE, ask=FALSE)"

RUN Rscript -e "remotes::install_version('alakazam', version='1.4.3', repos='https://cran.r-project.org'); \
                remotes::install_version('shazam',   version='1.3.2', repos='https://cran.r-project.org'); \
                remotes::install_version('scoper',   version='1.5.0', repos='https://cran.r-project.org'); \
                stopifnot(all(c('alakazam','shazam','scoper') %in% rownames(installed.packages())))"
