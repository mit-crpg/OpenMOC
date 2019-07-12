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

RUN apt-get install git swig python3-dev python3-numpy python3-matplotlib python3-h5py -y

#Build OpenMOC
WORKDIR "/OpenMOC"
RUN python3 setup.py install --user
RUN python3 setup.py install --user
