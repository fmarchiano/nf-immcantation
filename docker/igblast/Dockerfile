# pre-igblast: k-mer pre-filtered IgBLAST
#
# Build context: the pre-igblast source repo (not this directory).
#   docker build -t <tag> -f docker/igblast/Dockerfile /path/to/pre-igblast/

FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl gcc libc6-dev && \
    rm -rf /var/lib/apt/lists/*

RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install --no-cache-dir maturin

COPY rust/ /build/rust/
WORKDIR /build/rust
RUN maturin build --release --out /build/wheels

COPY pyproject.toml /build/
COPY src/ /build/src/
RUN pip wheel --no-deps --wheel-dir /build/wheels /build


FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl libc6 libgomp1 procps && \
    curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
        | tar xz -C /usr/local/bin s5cmd && \
    rm -rf /var/lib/apt/lists/*

RUN curl -fsSL https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz \
    | tar xz -C /opt && \
    ln -s /opt/ncbi-igblast-1.22.0/bin/igblastn /usr/local/bin/igblastn

COPY --from=builder /build/wheels /tmp/wheels
RUN pip install --no-cache-dir /tmp/wheels/*.whl awscli && rm -rf /tmp/wheels

ENV PATH="/opt/ncbi-igblast-1.22.0/bin:${PATH}"
