FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        procps curl \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir packaging presto==0.7.2 awscli
