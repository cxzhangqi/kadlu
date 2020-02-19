import numpy as np
import matplotlib
#matplotlib.use('svg')
matplotlib.use('Qt5Agg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import griddata
import cartopy
import cartopy.crs as ccrs
#import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import os
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources.chs import Chs
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.wwiii import Wwiii
import imageio


#terrain = cimgt.Stamen('terrain')
#terrain = cimgt.StamenTerrain()
#terrain = cimgt.Stamen('toner-lite')
#terrain = cimgt.MapQuestOpenAerial()

proj = ccrs.Mercator()
config = dict(
        bgcontour   = lambda v: np.linspace(min(v)-.1, max(v)+.1, 3),
        bathy=dict(
            cm      = plt.cm.bone,
            alpha   = 0.9,
            levels  = lambda v, n=25: np.linspace(min(v), max(v)-2, n),
            norm    = matplotlib.colors.Normalize(vmin=0, vmax=1000),
            title   = 'bathymetry (metres)'),
        temp=dict(
            cm      = plt.cm.coolwarm, 
            alpha   = 0.8,
            levels  = lambda v, n=12: np.linspace(min(v)+.1, max(v), n),
            norm    = matplotlib.colors.Normalize(vmin=-6, vmax=17),
            title   = 'temperature (celsius)'),
        salinity=dict(
            cm      = plt.cm.viridis,
            alpha   = 0.7,
            levels  = lambda v, n=12: np.linspace(min(v)+.1, max(v), n),
            norm    = matplotlib.colors.Normalize(vmin=20, vmax=40),
            title   = 'salinity (g/kg salt in water)'),
        waveheight=dict(
            #cm      = plt.cm.Spectral_r,
            #cm      = plt.cm.BrBG,
            cm      = plt.cm.BuPu,
            alpha   = 0.85,
            levels  = lambda v, n=12: np.linspace(min(v)+.1, max(v), n),
            norm    = matplotlib.colors.Normalize(vmin=0, vmax=15),
            title   = 'wave height (metres)')
    )



def plot2D(val, lat, lon, var, **kwargs): 
    """

    kwargs = dict(
            start=datetime(2015, 1, 9), end=datetime(2015, 1, 9, 3),
            south=44,                   west=-64.5, 
            north=45,                   east=-62.5, 
            top=0,                      bottom=0
        )

    kwargs = dict(
            start=datetime(2015, 1, 9), end=datetime(2015, 1, 10),
            south=45,                   west=-68.4, 
            north=51.5,                 east=-56.5, 
            top=0,                      bottom=5000
        )

    #Chs().fetch_bathymetry(**kwargs)
    #Hycom().fetch_temp(**kwargs)
    #Hycom().fetch_salinity(**kwargs)
    #Era5().fetch_windwaveswellheight(**kwargs)
    val, lat, lon, time = Era5().load_windwaveswellheight(**kwargs)
    var = 'waveheight'
    val, lat, lon, time, depth =  Hycom().load_temp(**kwargs)
    var = 'temp'
    val, lat, lon, time, depth =  Hycom().load_salinity(**kwargs)
    var = 'salinity'
    val, lat, lon =  Chs().load_bathymetry(**kwargs)
    var = 'bathy'

    """

    # project data onto coordinate space
    extent = proj.transform_points(
            ccrs.Geodetic(),
            np.array([kwargs['west'], kwargs['east']]), 
            np.array([kwargs['south'], kwargs['north']])
        )[:,:-1]
    projected_lonlat = proj.transform_points(
            ccrs.Geodetic(),
            lon,
            lat
        )
    plon = projected_lonlat[:,0]
    plat = projected_lonlat[:,1]
    # create interpolation grid, perform interpolation
    num_lats = 500
    num_lons = 500
    num_lats = 1000 
    num_lons = 1000 
    lons = np.linspace(start=min(plon), stop=max(plon), num=num_lons)
    lats = np.linspace(start=min(plat), stop=max(plat), num=num_lats)
    data = griddata(points=(plon, plat), values=val, xi=(lons[None,:],lats[:,None]), method='linear')
    coast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    fg = (.92, .92, .92, 1)
    fg = (0.66,0.66,0.66, 1)

    fname = f'{var}_{kwargs["start"].date().isoformat()}.png'
    print('saving ' + fname + '...')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, 
            title=config[var]['title']+f'\n{kwargs["start"].date().isoformat()}',
            projection=proj, 
            facecolor=config[var]['cm'](256), 
            frameon=True
        )
    ax.contourf(lons, lats, data,
                transform=proj,
                levels=config[var]['levels'](val),
                cmap=config[var]['cm'], 
                alpha=config[var]['alpha'],
                zorder=8
            )
    ax.contour(lons, lats, data,
                transform=proj,
                levels=config[var]['levels'](val),
                cmap=config[var]['cm'],
                alpha=1,
                linewidths=2,
                zorder=9
            )
    #ax.add_feature(ocean, facecolor=cm(0), edgecolor=(0,0,0,1), zorder=10)
    ax.add_feature(coast, facecolor=fg, edgecolor=(0,0,0,1), zorder=10)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--',
            zorder=11)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    ax.tick_params(axis='x', rotation=45)
    #ax.set_xticks(range(28, 35))
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=config[var]['norm'], 
                cmap=config[var]['cm']))

    plt.savefig(f'http/{var}/{fname}', bbox_inches='tight', 
                dpi=250, optimize=True)
    """
    plt.tight_layout()
    plt.show()
    """
    plt.close()

    return


def animate(kwargs, var, fetchfcn, loadfcn, debug=False):
    """

    var='temp'

    kwargs = dict(
            start=datetime(2015, 3, 1), end=datetime(2015, 9, 1),
            south=45,                   west=-68.4, 
            north=51.5,                 east=-56.5, 
            top=0,                      bottom=0
        )

    fetch = (Hycom().fetch_temp, Hycom().fetch_salinity, Era5().fetch_windwaveswellheight)
    load = (Hycom().load_temp, Hycom().load_salinity, Era5().load_windwaveswellheight)
    vartypes = ('temp', 'salinity', 'waveheight')
    vartypes = ('salinity', 'waveheight')


    for f,l,v in zip(fetch, load, vartypes):
        animate(kwargs, v, f, l)
    

    fetchfcn = Hycom().fetch_salinity
    loadfcn = Hycom().load_salinity

    """

    dirname = f'http/{var}'
    if not os.path.isdir(dirname): os.mkdir(dirname)
    png = lambda f: f if '.png' in f else None
    old = map(png, list(os.walk(f'http/{var}'))[0][2])
    _rm = [os.remove(f'http/{var}/{x}') for x in old if debug and x is not None]
    qry = kwargs.copy()
    t = datetime(qry['start'].year, 
                 qry['start'].month, 
                 qry['start'].day)


    while (t <= kwargs['end']):
        qry['start'] = t
        qry['end'] = t + timedelta(days=1)
        fname = f'http/{var}/{var}_{t.date().isoformat()}.png'
        if not os.path.isfile(fname):
            fetchfcn(**qry)
            if var in ('temp', 'salinity'): 
                val, lat, lon, time, depth = loadfcn(**qry)
            else:
                val, lat, lon, time = loadfcn(**qry)
            plot2D(val, lat, lon, var, **qry)
        t += timedelta(days=1)

    print(f'animating {var}.gif...')
    frames = sorted([f'http/{var}/{i}' for i in 
            map(png, list(os.walk(f'http/{var}'))[0][2]) 
            if  datetime.strptime(i, f'{var}_%Y-%m-%d.png') >= kwargs['start']
            and datetime.strptime(i, f'{var}_%Y-%m-%d.png') <  kwargs['end']])
    with imageio.get_writer(f'http/{var}.gif', mode='I') as w:
        list(map(w.append_data, map(imageio.imread, frames)))

    return



def fetch_topo():
    url = 'ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/shp/Elevation/canvec_50K_NS_Elevation_shp.zip'
    with requests.get(url, stream=True) as payload:
        assert payload.status_code == 200, 'error bad request'
        fname = storage_cfg()+'topo_NS.zip'
        with open(fname, 'wb') as f: f.write(payload.content)
        


#ax.add_feature(cartopy.feature.LAKES, edgecolor='black', projection=proj)
#ax.add_feature(cartopy.feature.RIVERS, edgecolor='black')
#ax.add_feature(rivers_lakes, facecolor='black', edgecolor='black', projection=proj)

#ax.set_xticks([-63.75, -63.5, -63.25, -63], crs=proj)
#ax.set_yticks([47.25, 47.5, 47.75, 48], crs=ccrs.Miller())

#ax.add_image(terrain, 9)
#ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'])
#ax.add_feature(rivers_lakes, facecolor=(0.89,0.92,0.94,1), edgecolor='black')
#ax.add_feature(coast, facecolor=(0.89,0.92,0.93,0.8), edgecolor=(0,0,0,0.6))

