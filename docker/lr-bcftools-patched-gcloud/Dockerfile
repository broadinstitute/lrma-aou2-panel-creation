FROM ubuntu:24.04
ARG BCFTOOLS_VERSION=1.23
ARG HTSLIB_VERSION=1.23
ENV BCFTOOLS_VERSION=$BCFTOOLS_VERSION
ENV HTSLIB_VERSION=$HTSLIB_VERSION
RUN apt-get update && \
    apt-get install -y build-essential autoconf automake make perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libperl-dev libgsl0-dev git curl
RUN mkdir /tmp/src && cd /tmp/src && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --branch $HTSLIB_VERSION && \
    git clone https://github.com/samtools/bcftools.git --branch $BCFTOOLS_VERSION
# Apply retry/resilience patch for libcurl HTTP backend (https://github.com/samtools/htslib/pull/1987)
RUN cd /tmp/src/htslib && \
    curl -L https://github.com/samtools/htslib/pull/1987.diff | git apply
# Apply silent truncation fix for synced reader I/O errors (https://github.com/samtools/bcftools/pull/2503)
RUN cd /tmp/src/bcftools && \
    curl -L https://github.com/samtools/bcftools/pull/2503.diff | git apply
RUN cd /tmp/src/htslib && autoreconf --install && ./configure --enable-gcs --enable-libcurl --enable-plugins && make -j8 && make install
RUN cd /tmp/src/bcftools && autoreconf --install && ./configure --enable-gcs --enable-configure-htslib && make -j8 && make install
RUN rm -rf /tmp/src

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 gnupg \
                 wget \
                 bc \
                 bedtools \
                 datamash \
                 gawk \
                 less \
                 tree \
                 vcftools && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-cli && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy purge gnupg && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    bcftools --version

ENV HOME=/root

# Define default command.
CMD ["bash"]
