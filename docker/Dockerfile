# syntax=docker/dockerfile:1

# Joseph S. Wirth
# May 2024

ARG VERSION=2.0.0

# start with ubuntu base
FROM ubuntu:22.04

# install python3, blast, wget, bowtie2, mummer, and cd-hit
RUN apt-get update && \
    apt-get -y install python3.11 python3-pip python3-setuptools python3-dev
RUN apt-get -y install ncbi-blast+
RUN apt-get -y install git
RUN apt-get -y install wget
RUN apt-get -y install mummer
RUN apt-get -y install bowtie2
RUN apt-get -y install cd-hit

# install sra-toolkit
RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar xzf sratoolkit.tar.gz -C /sra-toolkit
ENV PATH="${PATH}:/sra-toolkit"

# install ariba and biopython
RUN pip3 install Bio
RUN pip3 install ariba=2.12

# install espwAlleleCaller
RUN wget https://github.com/ncezid-biome/espwAlleleCaller/archive/refs/tags/v${VERSION}.tar.gz
RUN tar xzf v${VERSION}.tar.gz
RUN rm v${VERSION}.tar.gz
RUN mv /espwAlleleCaller-${VERSION} /espwAlleleCaller

# clone BIOME-scripts and copy miscDirectory.py file
RUN git clone https://github.com/ncezid-biome/BIOME-scripts.git
RUN echo "MISC_DIR = '/BIOME-scripts/misc-python-scripts/'" > /espwAlleleCaller/miscDirectory.py

# modify the executable so that it can be invoked directly from the command line
WORKDIR /espwAlleleCaller
RUN mv espwAlleleCaller.py espwAlleleCaller
RUN chmod +x espwAlleleCaller
ENV PATH="${PATH}:/espwAlleleCaller"

# make volume and set it as the workdir
VOLUME ["/data"]
WORKDIR /data

# set entry point to allow script execution outside the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
