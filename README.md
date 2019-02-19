# Welcome to Kadlu, a Python tookit for modelling underwater noise

Kadlu is under development, but will eventually 
contain a bunch of tools useful for modeling the underwater ocean 
soundscapes, for example:

 * Automated retrieval of relevant environmental data, including static 
   data such as bathymetry and seabed properties, and dynamic data such 
   as water temperature, salinity, and wave height.

 * Derivation of underwater acoustic properties (such as sound speed) from 
   the environmental data and conversion into format suitable for transmission 
   loss calculations.

 * Simulation of underwater noise produced by environmental forcings 
   such as waves and rain.

And potentially more ...

## Environmental noise model

Currently, we are working on translating [this transmission-loss code](https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu/tree/master/Nx2DSSFPE) 
from MATLAB to Python. As part of this work, we will be restructuring 
and documenting the code to make it more user-friendly and more easily 
adaptable to new scenarios.

## Dependencies and installation

Kadlu uses a number of standard Python libraries such as 
numpy, scipy, matplotlib, etc, which can be conveniently 
installed with pipe. Kadlu also uses a few C libraries:
 
  * [HDF5 (Hierarchical Data Format)](https://www.hdfgroup.org/) 
  * [NetCDF-4 (Network Common Data Form](https://www.unidata.ucar.edu/software/netcdf/)
  * [GDAL (Geospatial Data Abstraction Library)](https://www.gdal.org/)

which can be installed with apt-get.

The following bash script will install all the required 
libraries. It should work out of the box on newer Ubuntu systems.

```bash
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
```


## Notebook tutorials

 1. [Extract bathymetry data from a matlab file](docs/demo_notebooks/read_bathy.ipynb)

 2. [Polar and planar coordinates](docs/demo_notebooks/coordinates.ipynb)

 3. [Interpolate bathymetry data](docs/demo_notebooks/interp_bathy.ipynb)