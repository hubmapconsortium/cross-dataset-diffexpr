FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

## System packages
# I'm used to using `apt` but it complains that it doesn't have a stable
#   command-line interface, so let's use `apt-get` instead
RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    flex \
    git \
    libfreetype6-dev \
    libigraph0-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /root/.ssh

COPY bin /opt
COPY data /opt

WORKDIR /opt
