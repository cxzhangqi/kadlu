
import numpy as np
from kadlu.geospatial.data_sources import chs
from kadlu.geospatial.data_sources import wwiii
from kadlu.geospatial.data_sources import fetch_util
from datetime import datetime

storage_location = fetch_util.instantiate_storage_config()
filename = f"{storage_location}multi_1.glo_30m.hs.201702.grb2"
target_date = datetime(2017, 2, 3).replace(minute=0, hour=0, second=0, microsecond=0)

filename = f"{storage_location}multi_1.glo_30m.hs.201801.grb2"
target_date = datetime(2018, 1, 1).replace(minute=0, hour=0, second=0, microsecond=0)

######################

# mahone bay
south = 34.4
north = 44.7
west = -66.4
east = -61.8

# atlantic
north = 90
south = -90
west = -180
east = 180

north, south, east, west = None, None, None, None

grb = pygrib.open(filename)
#waveslice = grb.select(validDate=target_date, shortNameECMF='hs')[0]
waveslice = grb.select(name='Significant height of combined wind waves and swell')
#waveslice = grb.select(validDate=target_date)[0]

######################

lat = np.array([])
lon = np.array([])
data = np.array([])
for msg in grb:
    #mlat, mlon = np.array(msg.latlons())
    #vals, mlat, mlon = msg.data(0, False)
    vals, mlat, mlon = msg.data(south, north, west+180, east+180)
    vals[vals.mask] = 0
    #lat = np.append(lat, mlat[~vals.mask])
    #lon = np.append(lon, mlon[~vals.mask])
    #data = np.append(data, vals[~vals.mask])
    lat = np.append(lat, mlat)
    lon = np.append(lon, mlon)
    data = np.append(data, vals)
    break

grid = np.meshgrid(data, lat, lon, sparse=True)

######################

plt.figure()
m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='l')
#m=Basemap(projection='mill', lat_ts=10, llcrnrlon=west+180.0, urcrnrlon=east+179.5, llcrnrlat=lat.min(), urcrnrlat=lat.max(), resolution='c')
x, y = m(lon, lat)
#cs = m.pcolormesh(x, y, data, shading='flat', cmap=plt.cm.jet)
m.scatter(x, y, marker='.', color='r', zorder=3)
m.drawcoastlines(zorder=5)
m.fillcontinents()
m.drawmapboundary()
m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])
plt.colorbar(cs,orientation='vertical')
plt.tight_layout()
plt.title("testing")
plt.show()

