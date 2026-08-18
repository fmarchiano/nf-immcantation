# scoper-fast: Rust drop-in replacement for R SCOPer (spectral clustering)
#
# NOTE: This Dockerfile is a reference copy. To actually build the image,
# clone the scoper-fast source repo and build from there:
#
#   git clone <scoper-fast-repo-url>
#   cd scoper-fast
#   docker build -t fmarchiano/nf-immcantation-scoper-fast:0.1.0 -f docker/Dockerfile .
#
# The build context must be the scoper-fast repo root (not this directory),
# since it needs Cargo.toml, Cargo.lock, build.rs, src/, and bin/.

# ── Stage 1: builder ─────────────────────────────────────────────────
FROM rust:1.86-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        pkg-config \
        libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY src/ src/

RUN cargo build --release


# ── Stage 2: runtime ─────────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        libopenblas0 \
        ca-certificates \
        curl \
        procps \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/scoper-fast /usr/local/bin/
COPY bin/ /usr/local/bin/

CMD ["bash"]
