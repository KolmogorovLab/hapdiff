FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, fenderglass@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get -y install make gcc g++ && \
	apt-get -y install autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	#apt-get -y install samtools && \
	apt-get -y install python3-pip

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pysam scipy edlib matplotlib

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
	cd samtools-1.9/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

COPY . /opt/dipdiff
WORKDIR /opt/dipdiff
RUN make

ENV PATH "/opt/dipdiff:${PATH}"
ENV PYTHONUNBUFFERED "1"
ENV MPLCONFIGDIR "/tmp"
