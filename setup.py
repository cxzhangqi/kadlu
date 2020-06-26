from setuptools import setup, find_packages

import os

# create distribution and upload to pypi.org with:
#   $ python setup.py sdist bdist_wheel
#   $ twine upload dist/*

setup(name='kadlu',
        version=os.environ.get('KADLUVERSION', '0.0.0'), 
        description="MERIDIAN Python package for ocean ambient noise modelling",
        url='https://gitlab.meridian.cs.dal.ca/public_projects/kadlu',
        author='Oliver Kirsebom, Matthew Smith',
        author_email='oliver.kirsebom@dal.ca, matthew.smith@dal.ca',
        license='GNU General Public License v3.0',
        packages=find_packages(),
        install_requires=[
            'cartopy',
            #'cftime',
            'cdsapi',
            'geos',     # needed for cartopy
            'gsw',
            'imageio',
            'matplotlib',
            'mpl_scatter_density',
            'netcdf4',  # DEPENDS ON LIBRARIES:
            'numpy',
            'pandas',
            'Pillow',
            'proj',     # needed for cartopy
            'pygrib',   # DEPENDS ON eccodes 
            'pyproj',
            'pyqt5',
            'pytest',
            'scipy',
            #'tqdm',
            ],
        setup_requires=['pytest-runner',],
        tests_require=['pytest',],
        include_package_data=True,
        python_requires='>=3.8',
        #zip_safe=False
    )

