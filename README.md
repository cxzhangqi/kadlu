# Welcome to Kadlu, a Python tookit for modelling underwater noise

Kadlu is under development, but will eventually 
contain a bunch of tools useful for modeling the underwater ocean 
soundscapes, for example:

 * Automated retrieval of relevant environmental data, including static 
   data such as bathymetry and seabed properties, and dynamic data like 
   water temperature, salinity, and wave height.

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

## Dependencies

Kadlu uses a number of standard Python libraries such as 
numpy, scipy, and matplotlib, plus a few C libraries:
 
  * [HDF5 (Hierarchical Data Format)](https://www.hdfgroup.org/) 
  * [NetCDF-4 (Network Common Data Form)](https://www.unidata.ucar.edu/software/netcdf/)
  * [GDAL (Geospatial Data Abstraction Library)](https://www.gdal.org/)

Installation of these libraries is most easily accomplished using Anaconda.

## Installation with Anaconda
 
 1. [Download and install Anaconda](https://docs.anaconda.com/anaconda/install/).<br/>
    Make sure you get the Python 3.x version.<br/>
    Note that the Download button takes you to the macOS installer regardless of your OS, so make sure to pick the installer appropriate for your OS (Linux, macOS, Windows) 
 
 2. Clone the Kadlu repository
    ```terminal
      git clone https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu.git
      cd kadlu
    ```
    Note: Windows users may find it easier to download the source code using Gitlab's ''Download zip'' option.

 3. Create and activate Anaconda environment (this installs all dependencies)
    ```terminal
      conda env create -f environment.yml
      conda activate kadlu_env
    ```
 
 4. Install Kadlu
    ```terminal
      python setup.py sdist
      pip install dist/kadlu-0.0.1.tar.gz
    ```
 
 5. Check that everything is working by running pytest
    ```terminal
      pytest
    ```

## Installation without Anaconda

If you prefer to avoid using Anaconda, try to run this [bash script](https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu/blob/master/install_dep.sh). It will install all the required 
libraries using pip and apt-get. It should work out of the box on newer Ubuntu systems.


## Notebook tutorials

 1. [Extract bathymetry data from a matlab file](docs/demo_notebooks/read_bathy.ipynb)

 2. [Polar and planar coordinates](docs/demo_notebooks/coordinates.ipynb)

 3. [Interpolate bathymetry data](docs/demo_notebooks/interp_bathy.ipynb)

 4. [How to work with bathymetry data from the Canadian Hydrographic Service (GeoTIFF)](docs/demo_notebooks/CHS.ipynb)


## Useful resources

 *  [gsw Python package](https://github.com/TEOS-10/GSW-Python) (Python implementation of the Thermodynamic Equation of Seawater 2010)