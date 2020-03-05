Introduction
============

Kadlu is a software package for modeling underwater noise.
It was developed for the purpose of modeling noise due to waves and rain in shallow coastal waters, but contains tools useful for many other noise modeling tasks.

Kadlu is written in Python and utilizes a number of powerful software packages 
including NumPy, HDF5, NetCDF-4, SQLite, and GDAL.
It is licensed under the `GNU GPLv3 license <https://www.gnu.org/licenses/>`_ and hence freely available for anyone to use and modify.
The project is hosted on GitLab at 
`https://gitlab.meridian.cs.dal.ca/public_projects/kadlu <https://gitlab.meridian.cs.dal.ca/public_projects/kadlu>`_ .

Kadlu was developed by the `MERIDIAN <http://meridian.cs.dal.ca/>`_ Data Analytics Team at the 
`Institute for Big Data Analytics <https://bigdata.cs.dal.ca/>`_ at Dalhousie University with the 
support and assistance of David Barclay and Calder Robinson, both from the Department of Oceanography 
at Dalhousie University.

The first version of Kadlu, released in March 2020, provides functionalities for fetching environmental data (bathymetry, water temperature and salinity, wave height, wind speed, etc.) from online sources and loading into numpy arrays, interpolation on any coordinate array or grid, and plotting. 
Functionalities for sound propagation modelling will be included in the next release, which we anticipate will happen in May 2020.

The intended users of Kadlu are researchers and students in underwater acoustics working with ambient noise modeling. 
While Kadlu comes with complete documentation and comprehensive step-by-step tutorials, some familiarity with Python and 
especially the NumPy package would be beneficial. A basic understanding of 
the physical principles of underwater sound propagation would also be an advantage.

To get started with Kadlu, follow the <installation> and then proceed to the <tutorials/index> section.

In Inuit mythology, Kadlu is the name of a goddess that creates thundery weather, for example, by jumping on hollow ice. Thus, the name Kadlu was chosen to highlight the software package's main intended application, modeling of noise due to environmental forcings.

