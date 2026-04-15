FROM ubuntu:24.04
ARG BCFTOOLS_VERSION=1.23.1
ARG HTSLIB_VERSION=1.23.1
ENV BCFTOOLS_VERSION=$BCFTOOLS_VERSION
ENV HTSLIB_VERSION=$HTSLIB_VERSION
RUN apt-get update && \
    apt-get install -y build-essential autoconf automake make perl zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-openssl-dev libssl-dev libperl-dev libgsl0-dev git curl
RUN mkdir /tmp/src && cd /tmp/src && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --branch $HTSLIB_VERSION && \
    git clone https://github.com/samtools/bcftools.git --branch $BCFTOOLS_VERSION
# Apply retry/resilience patch for libcurl HTTP backend (https://github.com/samtools/htslib/pull/1987, still not merged into 1.23.1)
RUN cd /tmp/src/htslib && \
    curl -L https://github.com/samtools/htslib/pull/1987.diff | git apply
RUN cd /tmp/src/htslib && autoreconf --install && ./configure --enable-gcs --enable-libcurl --enable-plugins && make -j8 && make install
RUN cd /tmp/src/bcftools && autoreconf --install && ./configure --enable-gcs --enable-configure-htslib && make -j8 && make install
