FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, fenderglass@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get -y install make gcc g++ && \
	apt-get -y install samtools && \
	apt-get -y install python3-pip

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pysam scipy edlib matplotlib

COPY . /opt/dipdiff
WORKDIR /opt/dipdiff
RUN make

ENV PATH "/opt/dipdiff:${PATH}"
ENV PYTHONUNBUFFERED "1"
ENV MPLCONFIGDIR "/tmp"
