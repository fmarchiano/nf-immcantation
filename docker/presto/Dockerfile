FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir packaging presto==0.7.2
