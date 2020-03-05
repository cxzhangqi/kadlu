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
 
Check that everything is working by running pytest: ::

    pytest kadlu/ --doctest-modules
