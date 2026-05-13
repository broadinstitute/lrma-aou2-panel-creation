# python + bcftools + gcloud
# stage 0+1 from us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23.1 (https://github.com/broadinstitute/long-read-pipelines/blob/sh_update_gcloud_samtools_for_libcurl/docker/lr-gcloud-samtools/Dockerfile),
# but starting from python base, using gpg instead of deprecated apt-key, removing purges, and adding build-essential, gcc, and other -dev headers

############### stage 0: build samtools and bcftools from source
# FROM ubuntu:24.04 AS builder
FROM python:3.11-slim AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.23.1
ARG BCFTOOLS_VERSION=1.23.1
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 libdeflate-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} \
        && ./configure --without-curses --enable-libcurl --with-libdeflate \
        && make -j 4 -s all all-htslib \
        && make install install-htslib && \
    cd - && rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} \
        && ./configure --without-curses \
        && make -j 4 -s \
        && make install && \
    cd - && rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    # REMOVED: apt-get purge autoconf automake bzip2 gcc make wget
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help

############### stage 1: other commonly used bioinformatics utilities
# FROM ubuntu:24.04
# FIXED: Updated to standard Python to match Stage 0
FROM python:3.11-slim

ENV HOME=/root

# Define default command.
CMD ["bash"]

# copy from previous stage the binaries from samtools build
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/libexec/bcftools /usr/local/libexec/bcftools

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 build-essential \
                 gcc \
                 libcurl4t64 \
                 libdeflate-dev \
                 liblzma-dev \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 apt-transport-https \
                 ca-certificates \
                 gnupg \
                 curl \
                 wget \
                 bc \
                 bedtools \
                 datamash \
                 gawk \
                 less \
                 rename \
                 pigz \
                 tabix \
                 tree \
                 vcftools \
                 zlib1g-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-cli && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    # REMOVED: apt-get purge gnupg
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --version && \
    bcftools --version

COPY re-auth.sh /opt/
