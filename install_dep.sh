#!/bin/bash

# Update base container install
RUN apt-get update && apt-get upgrade -y

# Add unstable repo to allow us to access latest GDAL builds
echo deb http://ftp.uk.debian.org/debian unstable main contrib non-free >> /etc/apt/sources.list  
apt-get update

# Existing binutils causes a dependency conflict, correct version will be installed when GDAL gets intalled
apt-get remove -y binutils

# Install GDAL dependencies
apt-get -t unstable install -y libgdal-dev g++

# Update apt-get
apt-get install --reinstall build-essential && \
apt-get install software-properties-common && \
apt-get update

# Install HDF5 and NetCDF-4 C libraries with apt-get
apt-get install libhdf5-serial-dev libnetcdf-dev libnetcdff-dev

# Install tk with apt-get
apt-get install tk

# Update C env vars so compiler can find gdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal

# Install python packages specified in requirements.txt
pip install -r requirements.txt