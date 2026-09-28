FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        procps curl \
    && curl -fsSL https://github.com/peak/s5cmd/releases/download/v2.2.2/s5cmd_2.2.2_Linux-64bit.tar.gz \
       | tar xz -C /usr/local/bin s5cmd \
    && curl -fsSL https://github.com/torognes/vsearch/releases/download/v2.29.2/vsearch-2.29.2-linux-x86_64.tar.gz \
       | tar xz --strip-components=2 -C /usr/local/bin vsearch-2.29.2-linux-x86_64/bin/vsearch \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir packaging presto==0.7.2 awscli
