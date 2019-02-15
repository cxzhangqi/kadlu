# Use an official Python runtime as a parent image
FROM python:3.6-stretch

# Set the working directory to /app
WORKDIR /app

# Copy pyost package
COPY pyost/ /app/pyost
COPY requirements.txt /app
COPY setup.py /app

# Update base container install
RUN apt-get update && \
  apt-get upgrade -y

# Add unstable repo to allow us to access latest GDAL builds
RUN echo deb http://ftp.uk.debian.org/debian unstable main contrib non-free >> /etc/apt/sources.list  
RUN apt-get update

# Existing binutils causes a dependency conflict, correct version will be installed when GDAL gets intalled
RUN apt-get remove -y binutils

# Install GDAL dependencies
RUN apt-get -t unstable install -y libgdal-dev g++

# Install HDF5, netCDF4 and tk with apt-get.
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

# This will install latest version of GDAL
RUN pip install GDAL>=2.2.4

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

RUN python setup.py sdist
RUN pip install dist/pyost-0.0.1.tar.gz

# Run pytest when the container launches
ENTRYPOINT ["pytest"]

