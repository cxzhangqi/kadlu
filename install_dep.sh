#!/bin/bash

# Update base container install
sudo apt-get update
sudo apt-get install -y --reinstall build-essential
sudo apt-get install -y software-properties-common

# Install GDAL dependencies
sudo apt-get install -y gdal-bin libgdal-dev g++

# Install HDF5 and NetCDF-4 C libraries
sudo apt-get install -y libhdf5-serial-dev libnetcdf-dev libnetcdff-dev

# Install tk
sudo apt-get install -y tk

# Update C env vars so compiler can find gdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal

# Install python packages specified in requirements.txt
pip install -r requirements.txt