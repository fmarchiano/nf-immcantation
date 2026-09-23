# presto-fast: Rust-accelerated drop-in replacements for pRESTo tools
# Tools: MaskPrimersFast, AssemblePairsFast, CollapseSeqFast, PairSeqFast, BuildConsensusFast
#
# Build context: the presto source repo (not this directory).
#   docker build -t <tag> -f docker/presto-fast/Dockerfile /path/to/presto/

# ── Stage 1: builder ─────────────────────────────────────────────────
FROM python:3.11-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install --no-cache-dir maturin

WORKDIR /build

COPY maskprimers-fast/rust           maskprimers-fast/rust
COPY maskprimers-fast/src            maskprimers-fast/src
COPY maskprimers-fast/pyproject.toml maskprimers-fast/

COPY assemblepairs-fast/rust           assemblepairs-fast/rust
COPY assemblepairs-fast/src            assemblepairs-fast/src
COPY assemblepairs-fast/pyproject.toml assemblepairs-fast/

COPY collapseseq-fast/rust           collapseseq-fast/rust
COPY collapseseq-fast/src            collapseseq-fast/src
COPY collapseseq-fast/pyproject.toml collapseseq-fast/

COPY buildconsensus-fast/rust           buildconsensus-fast/rust
COPY buildconsensus-fast/src            buildconsensus-fast/src
COPY buildconsensus-fast/pyproject.toml buildconsensus-fast/

COPY pairseq-fast/rust           pairseq-fast/rust
COPY pairseq-fast/src            pairseq-fast/src
COPY pairseq-fast/pyproject.toml pairseq-fast/

# Build Rust extensions
RUN maturin build --release --manifest-path maskprimers-fast/rust/Cargo.toml --out /wheels
RUN maturin build --release --manifest-path assemblepairs-fast/rust/Cargo.toml --out /wheels
RUN maturin build --release --manifest-path collapseseq-fast/rust/Cargo.toml --out /wheels
RUN maturin build --release --manifest-path buildconsensus-fast/rust/Cargo.toml --out /wheels
RUN maturin build --release --manifest-path pairseq-fast/rust/Cargo.toml --out /wheels

# Install Rust wheels first, then Python packages (which depend on the Rust extensions)
RUN pip install --no-cache-dir /wheels/*.whl
RUN pip install --no-cache-dir \
        ./maskprimers-fast/ \
        ./assemblepairs-fast/ \
        ./collapseseq-fast/ \
        ./buildconsensus-fast/ \
        ./pairseq-fast/


# ── Stage 2: runtime ─────────────────────────────────────────────────
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps curl \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir awscli

COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin/MaskPrimersFast      /usr/local/bin/
COPY --from=builder /usr/local/bin/AssemblePairsFast    /usr/local/bin/
COPY --from=builder /usr/local/bin/CollapseSeqFast      /usr/local/bin/
COPY --from=builder /usr/local/bin/BuildConsensusFast.py /usr/local/bin/
COPY --from=builder /usr/local/bin/PairSeqFast.py       /usr/local/bin/
COPY --from=builder /usr/local/bin/SplitSeq.py          /usr/local/bin/

CMD ["bash"]
