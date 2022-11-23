FROM continuumio/miniconda3
MAINTAINER Fabio Cunial
ARG work_dir=/cromwell_root
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------
ARG gcloud_version=405.0.1
ARG htslib_version=1.16
ARG samtools_version=1.16.1
ARG bcftools_version=1.16
ARG minimap2_version=2.24
ARG k8_version=0.2.4
ARG hifiasm_version=0.16.1
# - PAV has releases, but for simplicity we just clone the repo to resolve a 
# dependency on SVPOP. Should be done more rigorously.
# ------------------------------------------------------------------------------


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        apt-transport-https \
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
        default-jdk \
        default-jre \
        git \
        python3-pip

# GSUTIL
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# MINIMAP2 + PAFTOOLS
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v${k8_version}/k8-${k8_version}.tar.bz2 | tar -jxf - \
    && cp k8-${k8_version}/k8-Linux /usr/local/bin/k8 \
    && rm -rf k8-${k8_version}/
RUN wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && rm minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && cd minimap2-${minimap2_version}_x64-linux \
    && cp minimap2 /usr/local/bin/ \
    && cp paftools.js /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf minimap2-${minimap2_version}_x64-linux \
    && minimap2 --help

# HIFIASM
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${hifiasm_version}.tar.gz \
    && tar -xzf ${hifiasm_version}.tar.gz \
    && rm -f ${hifiasm_version}.tar.gz \
    && cd hifiasm-${hifiasm_version} \
    && make \
    && cp hifiasm /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf hifiasm-${hifiasm_version} \
    && hifiasm -h

# PAV
COPY ./pav_environment.yml .
RUN conda env create -f ./pav_environment.yml
COPY pav_test_install.py .
SHELL ["conda", "run", "-n", "lr-pav", "/bin/bash", "-c"]
RUN echo "conda init bash" > ~/.bashrc
RUN python pav_test_install.py
RUN git clone --recursive https://github.com/EichlerLab/pav.git
# We fixed an empty CSV bug in rule <call_merge_haplotypes> of file
# <pav/rules/call.snakefile>.
COPY pav_call.snakefile ./pav/rules/call.snakefile
# We fixed an empty BED bug in rule <vcf_write_vcf> of file
# <pav/rules/vcf.snakefile>.
COPY pav_vcf.snakefile ./pav/rules/vcf.snakefile
