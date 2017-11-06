FROM ubuntu:zesty

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip /tmp/fastqc.zip
ADD https://github.com/lh3/bwa/archive/master.tar.gz /tmp/bwa.tar.gz
ADD https://github.com/samtools/htslib/archive/master.tar.gz /tmp/htslib.tar.gz
ADD https://github.com/samtools/samtools/archive/master.tar.gz /tmp/samtools.tar.gz
ADD https://github.com/samtools/bcftools/archive/master.tar.gz /tmp/bcftools.tar.gz
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD . /tmp/procread

RUN set -e \
      && ln -sf /bin/bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install autoconf default-jre gcc libbz2-dev liblzma-dev libncurses5-dev \
                            libz-dev make pbzip2 pigz python3.6 unzip \
      && apt-get clean

RUN set -e \
      && unzip -d /usr/local/src /tmp/fastqc.zip \
      && chmod +x /usr/local/src/FastQC/fastqc \
      && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc

RUN set -e \
      && tar xvf /tmp/bwa.tar.gz -C /usr/local/src \
      && mv /usr/local/src/bwa-master /usr/local/src/bwa \
      && cd /usr/local/src/bwa \
      && make \
      && ln -s /usr/local/src/bwa/bwa /usr/local/bin

RUN set -e \
      && tar xvf /tmp/htslib.tar.gz -C /usr/local/src \
      && mv /usr/local/src/htslib-master /usr/local/src/htslib \
      && cd /usr/local/src/htslib \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/samtools.tar.gz -C /usr/local/src \
      && mv /usr/local/src/samtools-master /usr/local/src/samtools \
      && cd /usr/local/src/samtools \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/bcftools.tar.gz -C /usr/local/src \
      && mv /usr/local/src/bcftools-master /usr/local/src/bcftools \
      && cd /usr/local/src/bcftools \
      && make \
      && make install

RUN set -e \
      && ln -s /usr/bin/python3.6 /usr/local/bin/python3 \
      && /usr/local/bin/python3 /tmp/get-pip.py \
      && pip install -U --no-cache-dir pip \
      && pip install -U --no-cache-dir /tmp/procread \
      && rm -rf /tmp/*

ENTRYPOINT ["procread"]
