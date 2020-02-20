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

 5. Generate an API token

    Kadlu can fetch environmental data from a variety of remote resources, including NOAA, ECCC and ECMWF. In order to access data from the ECMWF, it is necessary to first obtain and configure an API key, using the instructions here for [Windows](https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows), [Mac](https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+macOS) or [Linux](https://cds.climate.copernicus.eu/api-how-to).
 
 6. Check that everything is working by running pytest
    ```terminal
      pytest
    ```

## Installation without Anaconda

If you prefer to avoid using Anaconda, try to run this [bash script](https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu/blob/master/install_dep.sh). It will install all the required 
libraries using pip and apt-get. It should work out of the box on newer Ubuntu systems. Don't forget to [generate an API token](https://cds.climate.copernicus.eu/api-how-to) to gain access to ECMWF data after installing.

## Configuration
Kadlu allows configuration for how data is accessed and stored on your machine. These preferences are defined in kadlu/config.ini

 1. Data storage location

    By default, a storage folder will be created in Kadlu's root folder. A custom location can be configured by placing the following inside config.ini
    ```ini
      [storage]
      storage_location = /path/to/data/storage/
    ```

 2. ECMWF - CDS API Token

    Kadlu uses ECMWF's Era5 dataset as one of the optional data sources for wave height/direction/period and wind speed data.
    In order to access Era5 reanalysis data from the ECMWF, it is necessary to first obtain an API token.
    This can be obtained by registering an account and visiting [Copernicus API](https://cds.climate.copernicus.eu/api-how-to). Once logged in, your token will be displayed in the box under heading 'Install the CDS API key'.
    With your token, add the following to config.ini:
    ```ini
      [cdsapi]
      url = https://cds.climate.copernicus.eu/api/v2
      key = {YOUR_TOKEN_HERE}
    ```

## Notebook tutorials

 1. [Extract bathymetry data from a matlab file](docs/source/tutorials/read_bathy_tutorial/read_bathy_tutorial.ipynb)

 2. [Polar and planar coordinates](docs/source/tutorials/coordinates_tutorial/coordinates_tutorial.ipynb)

 3. [Interpolate bathymetry data](docs/source/tutorials/interp_bathy_tutorial/interp_bathy_tutorial.ipynb)

 4. [How to work with bathymetry data from the Canadian Hydrographic Service (GeoTIFF)](docs/source/tutorials/CHS_tutorial/CHS_tutorial.ipynb)
 
 5. [Calculate transmission loss](https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu/blob/master/docs/source/tutorials/calc_tl_tutorial/calc_tl_tutorial.ipynb)


## Useful resources

 *  [gsw Python package](https://github.com/TEOS-10/GSW-Python) (Python implementation of the Thermodynamic Equation of Seawater 2010)
