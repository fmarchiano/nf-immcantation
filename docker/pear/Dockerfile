FROM condaforge/mambaforge:24.3.0-0

LABEL maintainer="Fabio Marchiano <fabio.marchiano14@gmail.com>"
LABEL description="PEAR - Paired-End reAd mergeR v0.9.11"

RUN mamba install --yes \
        --channel bioconda \
        --channel conda-forge \
        pear=0.9.11 \
    && mamba clean --all --yes

CMD ["pear"]
