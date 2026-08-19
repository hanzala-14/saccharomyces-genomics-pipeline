FROM python:3.11-slim-bookworm

ARG SRA_VERSION=3.4.1
ENV PATH="/opt/sratoolkit/bin:${PATH}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        aria2 \
        ca-certificates \
        coreutils \
        procps \        
        curl \
        gzip \
        libbz2-1.0 \
        liblzma5 \
        libncurses6 \
        libxml2 \
        pigz \
        zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" \
       -o /tmp/sratoolkit.tar.gz \
    && tar -xzf /tmp/sratoolkit.tar.gz -C /opt \
    && mv "/opt/sratoolkit.${SRA_VERSION}-ubuntu64" /opt/sratoolkit \
    && rm /tmp/sratoolkit.tar.gz \
    && python --version \
    && aria2c --version | head -1 \
    && pigz --version \
    && fasterq-dump --version 2>&1 | head -1

WORKDIR /data
