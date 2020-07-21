import os
from setuptools import setup, find_packages

setup(name='kadlu',
        version=os.environ.get('KADLUVERSION', '0.0.0'), 
        description="MERIDIAN Python package for ocean ambient noise modelling",
        url='https://gitlab.meridian.cs.dal.ca/public_projects/kadlu',
        author='Oliver Kirsebom, Matthew Smith',
        author_email='oliver.kirsebom@dal.ca, matthew.smith@dal.ca',
        license='GNU General Public License v3.0',
        packages=find_packages(exclude=('tests',)),
        install_requires=[
            'cartopy',
            'cdsapi',
            'geos',     # needed for cartopy
            'gsw',
            'imageio',
            'matplotlib',
            'mpl_scatter_density',
            'netcdf4',  # DEPENDS ON binaries
            'numpy',
            'pandas',
            'Pillow',
            'proj',     # needed for cartopy
            'pygrib',   # DEPENDS ON eccodes binaries
            'pyproj',
            'pyqt5',
            'scipy',
            #'tqdm',
            ],
        setup_requires=['pytest-runner',],
        tests_require=['pytest','pytest-parallel'],
        include_package_data=True,
        exclude_package_data={'':['tests']},
        python_requires='>=3.8',
    )

