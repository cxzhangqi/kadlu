.. _installation_instructions:

Installation
=============

Kadlu is most easily installed using the Anaconda package manager.
Anaconda is freely available from `docs.anaconda.com/anaconda/install <https://docs.anaconda.com/anaconda/install/>`_. 
Make sure you get the Python 3.7 version and make sure to pick the installer appropriate for your OS (Linux, macOS, Windows) 

Clone the Kadlu repository: ::

    git clone https://gitlab.meridian.cs.dal.ca/public_projects/kadlu.git
    cd kadlu

Create and activate Anaconda environment: ::

    conda env create -f environment.yml
    conda activate kadlu_env
 
Install Kadlu: ::
    
    python setup.py sdist
    pip install dist/kadlu-1.0.0.tar.gz

Configuration: ::

1. Data storage location

By default, a storage folder will be created in Kadlu's root folder. Alternatively, a custom location can be configured by placing the following inside config.ini:

    [storage]
    storage_location = /path/to/data/storage/

2. ECMWF - CDS API Token

Kadlu uses ECMWF's Era5 dataset as one of the optional data sources for wave height/direction/period and wind speed data.
In order to access Era5 reanalysis data from the ECMWF, it is necessary to first obtain an API token.
This can be obtained by registering an account and visiting [Copernicus API](https://cds.climate.copernicus.eu/api-how-to). Once logged in, your token will be displayed in the box under heading 'Install the CDS API key'.
With your token, add the following to config.ini:

    [cdsapi]
    url = https://cds.climate.copernicus.eu/api/v2
    key = {YOUR_TOKEN_HERE}


Check that everything is working by running pytest: ::

    pytest kadlu/ --doctest-modules
