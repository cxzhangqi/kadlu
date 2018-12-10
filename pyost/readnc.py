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