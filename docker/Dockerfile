# syntax=docker/dockerfile:1

# Joseph S. Wirth
# September 2023

# start with ubuntu base
FROM ubuntu:22.04

# install python3, blast, wget, bowtie2, mummer, and cd-hit
RUN apt-get update && \
    apt-get -y install python3.11 python3-pip python3-setuptools python3-dev
RUN apt-get -y install ncbi-blast+
RUN apt-get -y install wget
RUN apt-get -y install mummer
RUN apt-get -y install bowtie2
RUN apt-get -y install cd-hit

# install ariba and bio python
RUN pip3 install Bio
RUN pip3 install ariba

# clone repository
RUN wget https://github.com/ncezid-biome/espwAlleleCaller/archive/refs/tags/v1.0.0.tar.gz
RUN tar xzf v1.0.0.tar.gz
RUN rm v1.0.0.tar.gz
RUN mv /espwAlleleCaller-1.0.0 /espwAlleleCaller

# modify the executable so that it can be invoked directly from the command line
WORKDIR /espwAlleleCaller
RUN echo "#!/usr/bin/python3" > espwAlleleCaller
RUN echo "" >> espwAlleleCaller
RUN cat espwAlleleCaller.py >> espwAlleleCaller
RUN chmod +x espwAlleleCaller
ENV PATH="${PATH}:/espwAlleleCaller"

# copy param file
COPY ./docker_params.py /espwAlleleCaller/auxillary/params.py

# make volume and set it as the workdir
VOLUME ["/data"]
WORKDIR /data

# set entry point to allow script execution outside the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
