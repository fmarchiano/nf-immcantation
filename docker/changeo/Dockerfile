FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget curl \
        libgomp1 \
        procps \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && rm -rf /var/lib/apt/lists/*

# IgBLAST 1.22.0 — only the igblastn binary is needed
RUN wget -q https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz \
    && tar -xzf ncbi-igblast-1.22.0-x64-linux.tar.gz \
    && mv ncbi-igblast-1.22.0/bin/igblastn /usr/local/bin/ \
    && rm -rf ncbi-igblast-1.22.0* \
    && rm -rf /var/lib/apt/lists/*

# External igblast databases mounted here at runtime
ENV IGDATA=/usr/local/share/igblast

RUN pip install --no-cache-dir changeo==1.3.0 awscli
