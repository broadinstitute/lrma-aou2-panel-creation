FROM continuumio/miniconda3:24.1.2-0
ARG work_dir=/shapeit5
WORKDIR ${work_dir}

RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        pkg-config \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python-dev \
        bsdmainutils \
        default-jre \
        default-jdk \
        unzip \
        python3-pycurl \
        python3-edlib 
        

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.19.1
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && autoheader \ 
    && autoconf \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# Boost (compile both static and shared libraries)
RUN wget https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2 \
    && tar --bzip2 -xf boost_1_73_0.tar.bz2 \
    && rm boost_1_73_0.tar.bz2 \
    && cd boost_1_73_0 \
    && ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=/usr/local/ \
    && ./b2 install --prefix=/usr/local/ -j16 link=shared,static \
    && cd ${work_dir} \
    && rm -rf boost_1_73_0 \
    && find /usr/local/lib -name "libboost*.so*" -o -name "libboost*.a" | head -20

RUN ldconfig -p | grep -E 'libz\.so|libbz2\.so|liblzma\.so|libm\.so|libpthread\.so|libcurl\.so' \
    && g++ --version \
    && gcc --version \
    && pkg-config --cflags --libs zlib liblzma libcurl || true \
    && find /usr/lib -name "libbz2*" -type f \
    && find /usr/include -name "bzlib.h" -type f

# Set library paths
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# Shapeit5
RUN git clone --recurse-submodules https://github.com/odelaneau/shapeit5.git \
    && cd shapeit5/phase_rare \
    && make \
        HTSSRC=/usr/local \
        HTSLIB_INC=/usr/local/include \
        HTSLIB_LIB=/usr/local/lib/libhts.a \
        BOOST_INC=/usr/local/include \
        BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a \
        BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a \
        CXXFLAGS="-O3 -L/usr/local/lib" \
        LDFLAGS="-L/usr/local/lib -lboost_serialization" \
        -j 16 \
    && ls -la ./bin/ \
    && if [ -f ./bin/phase_rare_static ]; then cp ./bin/phase_rare_static ${work_dir}/; elif [ -f ./bin/phase_rare ]; then cp ./bin/phase_rare ${work_dir}/; else echo "Error: phase_rare binary not found" && exit 1; fi \
    && cd ${work_dir} \
    && ls -la phase_rare* \
    && ./phase_rare --help

    
# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && bcftools --help


