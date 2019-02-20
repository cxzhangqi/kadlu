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

This [bash script](https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu/blob/master/install_dep.sh) will install all the required 
libraries. It should work out of the box on newer Ubuntu systems.


## Notebook tutorials

 1. [Extract bathymetry data from a matlab file](docs/demo_notebooks/read_bathy.ipynb)

 2. [Polar and planar coordinates](docs/demo_notebooks/coordinates.ipynb)

 3. [Interpolate bathymetry data](docs/demo_notebooks/interp_bathy.ipynb)

 4. [How to work with bathymetry data from the Canadian Hydrographic Service (GeoTIFF)](docs/demo_notebooks/CHS.ipynb)