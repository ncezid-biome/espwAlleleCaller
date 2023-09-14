# syntax=docker/dockerfile:1

# Joseph S. Wirth
# September 2023

# start with ubuntu base
FROM ubuntu:22.04

# install python3 and blast
RUN apt-get update
RUN apt-get -y install python3.11 python3-pip python3-setuptools python3-dev
RUN apt-get -y install ncbi-blast+

# install ariba and bio python
RUN pip3 install ariba
RUN pip3 install Bio

# clone repository
RUN git clone https://github.com/ncezid-biome/espwAlleleCaller.git

# modify the executable so that it can be invoked directly from the command line
WORKDIR /espwAlleleCaller
RUN cat "#!/usr/bin/python3" > espwAlleleCaller
RUN cat "" >> espwAlleleCaller
RUN cat espwAlleleCaller.py >> espwAlleleCaller
RUN chmod +x /espwAlleleCaller
ENV PATH="${PATH}:/espwAlleleCaller"

# copy param file
COPY ./docker_params.py /espwAlleleCaller/auxillary/params.py

# make volume and set it as the workdir
VOLUME ["/data"]
WORKDIR /data

# set entry point to allow script execution outside the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
