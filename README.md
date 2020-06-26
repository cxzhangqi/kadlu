# Welcome to Kadlu, a Python package for modelling ocean ambient noise

Kadlu was developed for the purpose of modelling noise due to waves and rain in shallow coastal 
waters, but contains tools useful for many other underwater sound modelling tasks.

Kadlu is written in Python and utilizes a number of powerful software packages 
including [NumPy](https://numpy.org/), [HDF5](https://www.hdfgroup.org/), 
[NetCDF-4](https://www.unidata.ucar.edu/software/netcdf/), and [SQLite](https://www.sqlite.org/index.html).
It is licensed under the [GNU GPLv3 license](https://www.gnu.org/licenses/) 
and hence freely available for anyone to use and modify.
The project is hosted on GitLab at 
[https://gitlab.meridian.cs.dal.ca/public_projects/kadlu](https://gitlab.meridian.cs.dal.ca/public_projects/kadlu). 
Kadlu was developed by the [MERIDIAN](http://meridian.cs.dal.ca/) Data Analytics Team at the 
[Institute for Big Data Analytics](https://bigdata.cs.dal.ca/) at Dalhousie University with the 
support and assistance of David Barclay and Calder Robinson, both from the Department of Oceanography 
at Dalhousie University.

Kadlu provides functionalities that automate the process of fetching and interpolating 
environmental data necessary to model ocean ambient noise levels (bathymetry, water temperature 
and salinity, wave height, wind speed, etc.). It also includes various routines that allow 
accurate estimates of noise source levels and transmission losses in realistic ocean environments.
You can find more information about the technical aspects of how sound propagation is modelled in 
Kadlu in [this note](docs/source/_static/kadlu_sound_propagation_note.pdf).

The intended users of Kadlu are researchers and students in underwater acoustics working with ambient noise modeling. 
While Kadlu comes with complete documentation and comprehensive step-by-step tutorials, some familiarity with Python and 
especially the NumPy package would be beneficial. A basic understanding of 
the physical principles of underwater sound propagation would also be an advantage.


## Installation

Kadlu is most easily installed using the Anaconda package manager.
Anaconda is freely available from [docs.anaconda.com/anaconda/install](https://docs.anaconda.com/anaconda/install/). 
Kadlu runs on the most recent stable version of Python 3. 

 1. Download dependency list and install dependencies using anaconda
    ```bash
    curl https://gitlab.meridian.cs.dal.ca/public_projects/kadlu/-/raw/master/environment.yml > environment.yml
    conda env create -f environment.yml python=3.8
    ```

 2. Activate the conda environment
    ```bash
    conda activate kadlu_env
    ```

 3. Install Kadlu
    ```bash
    pip install kadlu
    ```


## Configuration


#### Optionally set the storage directory

Kadlu allows configuration for where data is stored on your machine. By default, a folder 'kadlu_data' will be created in the home directory. To specify a custom location, run the following:

```python
import kadlu
kadlu.storage_cfg(setdir='/specify/desired/path/here/')
```


#### Optionally add an API token for fetching ERA5 data

Kadlu uses ECMWF's Era5 dataset as one of the data sources for wave height/direction/period and wind speed data.
By default, an API token is included with kadlu, but if you intend to make frequent use of the Era5 dataset, please consider obtaining your own token.
This can be obtained by registering an account at [Copernicus API](https://cds.climate.copernicus.eu/api-how-to). Once logged in, your token and URL will be displayed on the aforementioned webpage under heading 'Install the CDS API key'.
Additionally, you will need to accept the [Copernicus Terms of Use](https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products) to activate the token.

Configure Kadlu to use the token by executing:

```python
kadlu.era5_cfg(key="TOKEN_HERE", url="URL_HERE")
```


## Jupyter notebook tutorials

 1. [The Ocean Module](docs/source/tutorials/ocean_module_tutorial/ocean_module_tutorial.ipynb)

 2. [Fetch and Load Environmental Data](docs/source/tutorials/fetch_load_tutorial/fetch_load_tutorial.ipynb)

 3. [Interpolate Multi-Dimensional Data](docs/source/tutorials/interp_tutorial/interp_tutorial.ipynb)

 4. [Plot and Export Data](docs/source/tutorials/plot_export_tutorial/plot_export_tutorial.ipynb)

 5. [Transmission Loss](docs/source/tutorials/transm_loss_tutorial/transm_loss_tutorial.ipynb)


## Useful resources

 *  [gsw Python package](https://github.com/TEOS-10/GSW-Python) (Python implementation of the Thermodynamic Equation of Seawater 2010)

