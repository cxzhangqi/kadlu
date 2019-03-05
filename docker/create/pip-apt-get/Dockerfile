# Use an official Python runtime as a parent image
FROM python:3.6-stretch

# Set the working directory to /app
WORKDIR /app

# Copy kadlu package
COPY kadlu/ /app/kadlu
COPY requirements.txt /app
COPY setup.py /app

# Update base container install
RUN apt-get update

# Install GDAL dependencies
RUN apt-get install -y gdal-bin libgdal-dev g++

# Install HDF5, netCDF4 and tk
RUN apt-get install -y --reinstall build-essential && \
  apt-get install -y software-properties-common && \
  apt-get update && \
  apt-get install -y --no-install-recommends \
  libhdf5-serial-dev \
  libnetcdf-dev libnetcdff-dev \
  tk && \
  rm -rf /var/lib/apt/lists/*

# Update C env vars so compiler can find gdal
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

# Install packages specified in requirements.txt (including GDAL)
RUN pip install -r requirements.txt

RUN python setup.py sdist
RUN pip install dist/kadlu-0.0.1.tar.gz