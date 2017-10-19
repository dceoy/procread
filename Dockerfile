FROM python:latest

ADD https://github.com/lh3/bwa/archive/master.tar.gz /tmp/bwa.tar.gz
ADD https://github.com/samtools/htslib/archive/master.tar.gz /tmp/htslib.tar.gz
ADD https://github.com/samtools/samtools/archive/master.tar.gz /tmp/samtools.tar.gz
ADD . /tmp/procread

RUN set -e \
      && ln -sf /bin/bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install pigz \
      && apt-get clean

RUN set -e \
      && tar xvf /tmp/bwa.tar.gz -C /usr/local/src \
      && cd /usr/local/src/bwa-master \
      && make \
      && ln -s /usr/local/src/bwa-master/bwa /usr/local/bin

RUN set -e \
      && tar xvf /tmp/htslib.tar.gz -C /usr/local/src \
      && cd /usr/local/src/htslib-master \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/samtools.tar.gz -C /usr/local/src \
      && cd /usr/local/src/samtools-master \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && pip install -U pip cutadapt \
      && pip install -U /tmp/procread \
      && rm -rf /tmp/*

ENTRYPOINT ["procread"]
