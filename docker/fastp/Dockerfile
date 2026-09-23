FROM condaforge/mambaforge:24.3.0-0

LABEL maintainer="Fabio Marchiano <fabio.marchiano14@gmail.com>"
LABEL description="fastp v1.3.3"

RUN mamba install --yes \
        --channel bioconda \
        --channel conda-forge \
        fastp=1.3.3 \
    && mamba clean --all --yes

RUN mamba install --yes --channel conda-forge awscli && mamba clean --all --yes

RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && apt-get purge -y curl && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

CMD ["fastp"]
