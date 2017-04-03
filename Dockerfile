FROM ubuntu:16.04

ENV DEBIAN_FRONTEND noninteractive

RUN mkdir /app

COPY . /app/

RUN apt-get update -y

RUN apt-get install software-properties-common -y

RUN add-apt-repository -y ppa:staticfloat/juliareleases
RUN add-apt-repository -y ppa:staticfloat/julia-deps

RUN apt-get update && apt-get upgrade -y


RUN apt-get install build-essential julia gfortran pkg-config wget libnlopt0 unzip -y

ENV MAKEFLAGS -j4

RUN cd /app && julia bin/INSTALL.jl;

WORKDIR /app

CMD ["julia bin/runInference.jl --help"]