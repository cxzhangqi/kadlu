from kadlu.bathy_reader import BathyReader, LatLon, write_bathy
from kadlu.bathy_interpolator import BathyInterpolator
import matplotlib.pyplot as plt

# folder where the GeoTIFF files are stored
folder = '../kadlu/tests/assets/tif/'

# reader
reader = BathyReader(folder) 

# interpolator (select between 'linear', 'nearest', 'cubic')
interp = BathyInterpolator(reader, method='nearest', latlon_SW=LatLon(42.1,-60.2), latlon_NE=LatLon(44.4,-59.1)) 

# evaluate at 43.2N 59.8W
depth = interp.eval_ll(lat=43.2, lon=-59.8)
print('Depth at 43.2N59.8W:  {0:.3f} km'.format(-depth/1e3))

# evaluate on a 2x2 grid and save to a MATLAB file
lats = [43.2, 43.4]
lons = [-59.8, -59.3]
depths = interp.eval_ll(lat=lats, lon=lons, grid=True)
write_bathy(lat=lats, lon=lons, bathy=depths, destination='bathy.mat')

# plot the bathymetry using a 200 x 200 grid
interp.plot_ll(lat_bins=200, lon_bins=200)

# draw the plot
plt.show()
