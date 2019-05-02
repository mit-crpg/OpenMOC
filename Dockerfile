#Use latest Ubuntu version as parent image
# FIXME: Decide on a specific version of ubuntu
FROM ubuntu:latest
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git && \
    apt-get install --reinstall build-essential -y

#Clone Repository
RUN git clone -b 3D-MOC https://github.com/mit-crpg/OpenMOC.git /OpenMOC/

#Install necessary python dependencies
RUN apt-get install python3.7 -y

#Response to Tzdata package
ENV TZ=US
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get install git swig python-dev python-numpy python-matplotlib python-h5py -y

#Build OpenMOC
WORKDIR "/OpenMOC"
RUN python setup.py install
