# Welcome to Kadlu, a Python package for modelling ocean ambient noise

Kadlu was developed for the purpose of modelling noise due to waves and rain in shallow coastal 
waters, but contains tools useful for many other underwater sound modelling tasks.

Kadlu is written in Python and utilizes a number of powerful software packages 
including [NumPy](https://numpy.org/), [HDF5](https://www.hdfgroup.org/), 
[NetCDF-4](https://www.unidata.ucar.edu/software/netcdf/), 
[SQLite](https://www.sqlite.org/index.html), and [GDAL](https://www.gdal.org/).
It is licensed under the [GNU GPLv3 license](https://www.gnu.org/licenses/) 
and hence freely available for anyone to use and modify.
The project is hosted on GitLab at 
[https://gitlab.meridian.cs.dal.ca/public_projects/kadlu](https://gitlab.meridian.cs.dal.ca/public_projects/kadlu). 
Kadlu was developed by the [MERIDIAN](http://meridian.cs.dal.ca/) Data Analytics Team at the 
[Institute for Big Data Analytics](https://bigdata.cs.dal.ca/) at Dalhousie University with the 
support and assistance of David Barclay and Calder Robinson, both from the Department of Oceanography 
at Dalhousie University.

You can install the latest version of Kadlu from the Python Package Index (PyPI) repository using Anaconda. 
For more information, please consult [Kadlu's Documentation Page](https://docs.meridian.cs.dal.ca/kadlu/).

The first version of Kadlu, released in March 2020, provides functionalities for fetching environmental data 
(bathymetry, water temperature and salinity, wave height, wind speed, etc.) from online sources and loading into 
numpy arrays, interpolation on any coordinate array or grid, and plotting. Functionalities for sound propagation 
modelling will be included in the next release, anticipated for May 2020.

The intended users of Kadlu are researchers and students in underwater acoustics working with ambient noise modeling. 
While Kadlu comes with complete documentation and comprehensive step-by-step tutorials, some familiarity with Python and 
especially the NumPy package would be beneficial. A basic understanding of 
the physical principles of underwater sound propagation would also be an advantage.


## Installation

Kadlu is most easily installed using the Anaconda package manager.
Anaconda is freely available from [docs.anaconda.com/anaconda/install](https://docs.anaconda.com/anaconda/install/). 
Make sure you get the Python 3.7 version and make sure to pick the installer appropriate for your OS (Linux, macOS, Windows) 

 1. Clone the Kadlu repository:
    ```bash
    git clone https://gitlab.meridian.cs.dal.ca/public_projects/kadlu.git
    cd kadlu
    ```

 2. Create and activate Anaconda environment:
    ```bash
    conda env create -f environment.yml
    conda activate kadlu_env
    ```

 3. Install the PyPI package manager and Jupyter Notebook:
    ```bash
    conda install pip
    conda install jupyter
    ```

 4. Configure Kadlu (see below)

 5. Install Kadlu:
    ```bash
    python setup.py sdist
    pip install dist/kadlu-1.0.0.tar.gz
    ```


## Configuration

Kadlu allows configuration for how data is accessed and stored on your machine. These preferences are defined in kadlu/config.ini

 1. Data storage location

    By default, a folder 'kadlu_data' will be created in the user's home directory. To specify a custom location, run the following code::
    ```python
    from kadlu import data_util
    data_util().storage_cfg(setdir='/specify/desired/path/here/')
    ```

 2. ECMWF - CDS API Token

    Kadlu uses ECMWF's Era5 dataset as one of the optional data sources for wave height/direction/period and wind speed data.
    In order to access Era5 reanalysis data from the ECMWF, it is necessary to first obtain an API token.
    This can be obtained by registering an account at [Copernicus API](https://cds.climate.copernicus.eu/api-how-to). Once logged in, your token and URL will be displayed on the aforementioned webpage under heading 'Install the CDS API key'.
    Configure Kadlu to use the token by executing:
    ```python
    from kadlu import data_util
    data_util().era5_cfg(key="TOKEN_HERE", url="URL_HERE")
    ```


## Check your installation

Check that everything is working by running pytest:
```bash
pytest kadlu/ --doctest-modules
```


## Jupyter notebook tutorials

 1. [The Ocean Module](docs/source/tutorials/ocean_module_tutorial/ocean_module_tutorial.ipynb)

 2. [Fetch and Load Environmental Data](docs/source/tutorials/fetch_load_tutorial/fetch_load_tutorial.ipynb)

 3. [Interpolate Multi-Dimensional Data](docs/source/tutorials/interp_tutorial/interp_tutorial.ipynb)

 4. [Plot and Export Data](docs/source/tutorials/plot_export_tutorial/plot_export_tutorial.ipynb)


## Useful resources

 *  [gsw Python package](https://github.com/TEOS-10/GSW-Python) (Python implementation of the Thermodynamic Equation of Seawater 2010)
