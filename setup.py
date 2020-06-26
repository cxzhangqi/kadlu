from setuptools import setup, find_packages

import os

# create distribution and upload to pypi.org with:
#   $ python setup.py sdist bdist_wheel
#   $ twine upload dist/*

setup(name='kadlu',
        #version='2.1.1',
        version=os.environ.get('KADLUVERSION', '0.0.0'), 
        description="MERIDIAN Python package for ocean ambient noise modelling",
        url='https://gitlab.meridian.cs.dal.ca/public_projects/kadlu',
        author='Oliver Kirsebom, Matthew Smith',
        author_email='oliver.kirsebom@dal.ca, matthew.smith@dal.ca',
        license='GNU General Public License v3.0',
        packages=find_packages(),
        install_requires=[
            'numpy',
            'scipy',
            'pytest',
            'cdsapi',
            'pyproj',
            'matplotlib',
            'geos',     # needed for cartopy
            'proj',     # needed for cartopy
            'Pillow',
            'imageio',
            'netcdf4',
            'mpl_scatter_density',
            'pygrib',
            'cartopy',
            'pyqt5',
            ],
        setup_requires=['pytest-runner',],
        tests_require=['pytest',],
        include_package_data=True,
        python_requires='>=3.8',
        #zip_safe=False
    )

