Installation
=============

Kadlu is most easily installed using the Anaconda package manager.
Anaconda is freely available from: ::

    https://docs.anaconda.com/anaconda/install/
    
Make sure you get the Python 3.x version. Also, note that the Download button takes you to the macOS installer regardless of your OS, so make sure to pick the installer appropriate for your OS (Linux, macOS, Windows) 

Clone the Kadlu repository: ::

    git clone https://gitlab.meridian.cs.dal.ca/data_analytics_dal/packages/kadlu.git
    cd kadlu

Create and activate Anaconda environment: ::

    conda env create -f environment.yml
    conda activate kadlu_env
 
Install Kadlu: ::
    
    python setup.py sdist
    pip install dist/kadlu-0.0.1.tar.gz
 
Check that everything is working by running pytest: ::

    pytest
