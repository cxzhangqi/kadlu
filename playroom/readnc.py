
# The netcdf4-python package is used to read 
# netCDF files. The packages can be obtained 
# from:
#
#    http://unidata.github.io/netcdf4-python/
#
# Note that netcdf4-python has a number of 
# requirements including the C libraries 
# HDF5 and netCDF-4. 

from netCDF4 import Dataset

d = Dataset('../res/GEBCO_2014_2D.nc')

print(d.file_format)

print(d.dimensions.keys())
print(d.variables.keys())

print(d.dimensions['lat'])
print(d.dimensions['lon'])

b = d.variables['elevation']
print(b.shape)
print(b[0:10])

import matplotlib.pyplot as plt
plt.plot(d.variables['elevation'][int(1.5*10800),:])
plt.show()
