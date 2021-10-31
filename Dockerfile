FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, fenderglass@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get -y install make gcc g++ && \
	apt-get -y install samtools && \
	apt-get -y install python3-pip

#3RUN apt-get update && \
#       apt-get -y install time git make wget autoconf gcc g++ && \
##      apt-get -y install autoconf bzip2 lzma-dev zlib1g-dev && \
#       apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev && \
#       apt-get -y install liblzma-dev libhdf5-dev libncurses5-dev libncursesw5-dev && \
#       apt-get -y install python3-dev python3-pip && \
#       apt-get -y install cmake && \
#       apt-get clean && \
#       apt-get purge && \
#       rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pysam scipy edlib matplotlib

COPY . /opt/dipdiff
WORKDIR /opt/dipdiff
RUN make

ENV PATH "/opt/dipdiff:${PATH}"
ENV PYTHONUNBUFFERED "1"
ENV MPLCONFIGDIR "/tmp"
