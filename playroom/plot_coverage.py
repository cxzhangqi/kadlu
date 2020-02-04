import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
matplotlib.use("TkAgg")
from scipy.interpolate import griddata
from kadlu.geospatial.data_sources.data_util import index
from kadlu.geospatial.data_sources.chs import Chs
from datetime import datetime

###

kwargs = dict(
        start=datetime(2015, 1, 9), end=datetime(2015, 1, 10, 12),
        south=44,                   west=-64.5, 
        north=46,                   east=-62.5, 
        top=0,                      bottom=5000
    )

Chs().fetch_bathymetry(**kwargs)

val1, lat1, lon1 =  np.array(Chs().load_bathymetry(**kwargs)).astype(float)

###

#plt.ion()
#plt.ioff()
interp=True

mapargs = dict(
        projection='mill',
        #projection='cyl',
        #projection='cea',  # cylindrical equal area
        lat_ts=10,
        #lat_ts=0,
        llcrnrlon=kwargs['west'],
        urcrnrlon=kwargs['east'],
        llcrnrlat=kwargs['south'],
        urcrnrlat=kwargs['north'],
        resolution='f'
    )

num_lats = int(np.ceil((kwargs['north']- kwargs['south']) / 0.001)) + 1
num_lons = abs(int(np.ceil((kwargs['west']- kwargs['east']) / 0.001)) + 1)
xi = np.linspace(start=min(lon1), stop=max(lon1), num=num_lons)
yi = np.linspace(start=min(lat1), stop=max(lat1), num=num_lats)
zi = griddata(points=(lon1, lat1), values=val1, xi=(xi[None,:],yi[:,None]), method='linear')
X, Y = np.meshgrid(xi, yi)
x, y = m(X, Y)
serialize(kwargs, (x, y, zi), 'plot_interp')
serialize(mapargs, Basemap(**mapargs))

plt.style.available
plt.style.use('grayscale')
plt.style.use('seaborn-deep')

fig = plt.figure()
m = deserialize(mapargs)
m.drawcoastlines()
m.fillcontinents()
m.drawmapboundary()
m.drawrivers()
m.drawparallels(np.arange(kwargs['south'], kwargs['north'], .25), labels=[1,0,0,0])
m.drawmeridians(np.arange(kwargs['west'],  kwargs['east'],  .25), labels=[0,0,0,1])
if interp:
    x, y, zi = deserialize(kwargs, 'plot_interp')
    plt.contourf(x, y, zi, cmap='Greys_r')
    plt.colorbar()
else:
    x, y = m(lon1, lat1)
    color = (val1 / -261)
    color[color > 1] = 1
    plt.scatter(x, y, c=color, s=30, marker='.', cmap='Greys_r')#, zorder=10)
    norm = matplotlib.colors.Normalize(vmin=min(val1), vmax=max(val1))
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap='Greys_r'))
fig.tight_layout(w_pad=0.05)
plt.title('depth (metres)')
plt.show()

