.. _installation_instructions:

Installation
=============

Kadlu is most easily installed using the Anaconda package manager.
Anaconda is freely available from `docs.anaconda.com/anaconda/install <https://docs.anaconda.com/anaconda/install/>`_. 
Make sure you get the Python 3 version and make sure to pick the installer appropriate for your OS (Linux, macOS, Windows) 

Update your Anaconda installation to the latest version of python: ::

    conda install -c anaconda python=3.8

Clone the Kadlu repository: ::

    git clone https://gitlab.meridian.cs.dal.ca/public_projects/kadlu.git
    cd kadlu

Create and activate Anaconda environment: ::

    conda env create -f environment.yml
    conda activate kadlu_env
 
Install the PyPI package manager and Jupyter Notebook: ::
    
    conda install pip
    conda install jupyter

Configure Kadlu: ::

0. Import Kadlu

    import kadlu

1. Data storage location

By default, a folder 'kadlu_data' will be created in the user's home directory. To specify a custom location, run the following code:

    kadlu.storage_cfg(setdir='/specify/desired/path/here/')

2. ECMWF - CDS API Token

Kadlu uses ECMWF's Era5 dataset as one of the optional data sources for wave height/direction/period and wind speed data.
In order to access Era5 reanalysis data from the ECMWF, it is necessary to first obtain an API token.
This can be obtained by registering an account at [Copernicus API](https://cds.climate.copernicus.eu/api-how-to). Once logged in, your token and URL will be displayed on the aforementioned webpage under heading 'Install the CDS API key'.
Additionally, you will need to accept the [Copernicus Terms of Use](https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products) to activate the token.
Configure Kadlu to use the token by executing:

    kadlu.era5_cfg(key="TOKEN_HERE", url="URL_HERE")

Install Kadlu: ::
    
    python setup.py sdist
    pip install dist/kadlu-2.0.0.tar.gz

Check that everything is working by running pytest: ::

    pytest kadlu/ --doctest-modules
