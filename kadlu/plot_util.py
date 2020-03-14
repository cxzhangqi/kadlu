import os
from datetime import datetime, timedelta
from multiprocessing import Process, Queue

import numpy as np
import imageio
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
#matplotlib.use('Qt5Agg')
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata

from kadlu.geospatial.data_sources.chs import Chs
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.wwiii import Wwiii
from kadlu.geospatial.data_sources.source_map import load_map
from kadlu.geospatial.data_sources.data_util import storage_cfg
from kadlu.geospatial.data_sources.fetch_handler import fetch_handler


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


def plot2D(var, source, plot_wind=False, save=False, **kwargs): 
    """

    kwargs = dict(
            start=datetime(2015, 8, 9), end=datetime(2015, 8, 10),
            south=45,                   west=-68.4, 
            north=51.5,                 east=-56.5, 
            top=0,                      bottom=10
        )

    plot2D('temp', 'hycom', save=True, **kwargs)

    #Chs().fetch_bathymetry(**kwargs)
    #Hycom().fetch_temp(**kwargs)
    Hycom().fetch_salinity(**kwargs)
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

    if f'{var}_{source}' not in load_map.keys():
        raise KeyError(f'could not find source for variable. valid vars and '
                       f'sources: {[k.split("_") for k in load_map.keys()]}')

    if 'start' not in kwargs.keys():
        kwargs['start'], kwargs['end'] = datetime.now(), datetime.now()

    loadfcn = load_map[f'{var}_{source}']
    data = loadfcn(**kwargs)
    val, lat, lon = data[:3]

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
    num_lats = 1000 
    num_lons = 1000 
    lons = np.linspace(start=min(plon), stop=max(plon), num=num_lons)
    lats = np.linspace(start=min(plat), stop=max(plat), num=num_lats)
    data = griddata(points=(plon, plat), values=val, xi=(lons[None,:],lats[:,None]), method='linear')
    coast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    fg = (.92, .92, .92, 1)
    #fg = (0.66,0.66,0.66, 1)
    fname = f'{var}_{kwargs["start"].date().isoformat()}.png'
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

    if plot_wind is not False:
        if plot_wind.lower() == 'era5': 
            windfcnU, windfcnV = (Era5().load_wind_u, Era5().load_wind_v)
        elif plot_wind.lower() == 'wwiii': 
            windfcnU, windfcnV = (Wwiii().load_wind_u, Wwiii().load_wind_v)
        else: 
            raise ValueError('invalid wind source. must be \'era5\' or \'wwiii\'')

        uval, ulat, ulon, utime = windfcnU(**kwargs)
        vval, vlat, vlon, vtime = windfcnV(**kwargs)
        assert(len(vval) == len(uval))  # this can be fixed with an SQL JOIN in load module
        ax.quiver(ulon, ulat, uval, vval, transform=ccrs.PlateCarree(), 
                regrid_shape=20, zorder=10)

    ax.add_feature(coast, facecolor=fg, edgecolor=(0,0,0,1), zorder=11)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--',
            zorder=11)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    ax.tick_params(axis='x', rotation=45)
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=config[var]['norm'], 
                cmap=config[var]['cm']))

    if save is not False:
        if not os.path.isdir(f'{storage_cfg()}figures'): 
            os.mkdir(f'{storage_cfg()}figures')
        print(f'saving figure to {storage_cfg()}figures/{fname if save is True else save}')
        plt.savefig(f'{storage_cfg()}figures/{fname if save is True else save}', 
                bbox_inches='tight', dpi=200, figsize=(12,8), optimize=True)
        plt.close()
    else: 
        plt.show()

    return


def animate(var, source, kwargs, step=timedelta(hours=12), fps=30, plot_wind=False, debug=False):
    """
    animate('temp', 'hycom', kwargs, step=timedelta(hours=6), fps=60, debug=False)

    var='temp'

    kwargs = dict(
            start=datetime(2015, 1, 2, 12), end=datetime(2015, 12, 31),
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
    
    """
    # download all the data
    fetch_handler(var, source, **kwargs)

    # prepare folder and check for existing frames
    dirname = storage_cfg() + 'figures/'
    if not os.path.isdir(dirname): os.mkdir(dirname)
    png = lambda f: f if '.png' in f else None
    old = map(png, list(os.walk(dirname))[0][2])
    _rm = [os.remove(f'{dirname}{x}') for x in old if debug and x is not None]

    # generate image frames
    qry = kwargs.copy()
    cur = datetime(kwargs['start'].year, kwargs['start'].month, kwargs['start'].day)
    while cur <= kwargs['end']:
        qry['start'] = cur
        qry['end'] = cur + step
        fname = f'{var}_{cur.isoformat()}.png'
        if not os.path.isfile(f'{dirname}/{fname}'): 
            plot2D(var, source, plot_wind=plot_wind, save=fname, **qry)
        cur += step

    print(f'animating {var}.gif...')
    fmt = f'{var}_%Y-%m-%dT%H:%M:%S.png'
    frames = sorted([f'{dirname}{i}' for i in 
            map(png, list(os.walk(f'{dirname}'))[0][2]) if i is not None
            and datetime.strptime(i, fmt) >= kwargs['start']
            and datetime.strptime(i, fmt) <= kwargs['end']])

    #with imageio.get_writer(f'http/{var}.gif', mode='I') as w:
    savedir = f'{storage_cfg()}animated'
    if not os.path.isdir(savedir): os.mkdir(savedir)

    with imageio.get_writer(
            f'{savedir}/{var}_{kwargs["start"].date().isoformat()}_{kwargs["end"].date().isoformat()}.mp4', 
            mode='I', macro_block_size=4, format='FFMPEG', fps=fps) as w:
        list(map(w.append_data, map(imageio.imread, frames)))

    return 


    """
https://imageio.readthedocs.io/en/stable/examples.html#writing-videos-with-ffmpeg-and-vaapi

    w = imageio.get_writer(f'{dirname}/salinity.mp4', 
            format='FFMPEG', mode='I', fps=30,
            codec='h264_vaapi',
            output_params=['-vaapi_device',
                           '/dev/dri/renderD128',
                           '-vf',
                           'format=gray|nv12,hwupload'],
            pixelformat='vaapi_vld')

    with imageio.get_writer(f'{dirname}/{var}.mp4', mode='I', format='FFMPEG') as w:
        list(map(w.append_data, map(imageio.imread, frames)))

    for f in frames[250:]:
        fdat = imageio.imread(f)
        print(fdat.shape, f)

    for f in frames[256:]:
        im = Image.open(f)
        im = im.resize((1085, 839))
        im.save(f)
        im.close()

    """

"""
from ridge_map import RidgeMap
RidgeMap().plot_map()
rmap = RidgeMap((kwargs['west'], kwargs['south'], kwargs['east'], kwargs['north']))
rvals = rmap.get_elevation_data(num_lines=150)
rmap.plot_map(label='', ax=ax)
plt.show()
"""


"""
def fetch_topo():
    url = 'ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/shp/Elevation/canvec_50K_NS_Elevation_shp.zip'
    with requests.get(url, stream=True) as payload:
        assert payload.status_code == 200, 'error bad request'
        fname = storage_cfg()+'topo_NS.zip'
        with open(fname, 'wb') as f: f.write(payload.content)
"""
        


#terrain = cimgt.Stamen('terrain')
#terrain = cimgt.StamenTerrain()
#terrain = cimgt.Stamen('toner-lite')
#terrain = cimgt.MapQuestOpenAerial()


#ax.add_feature(cartopy.feature.LAKES, edgecolor='black', projection=proj)
#ax.add_feature(cartopy.feature.RIVERS, edgecolor='black')
#ax.add_feature(rivers_lakes, facecolor='black', edgecolor='black', projection=proj)

#ax.set_xticks([-63.75, -63.5, -63.25, -63], crs=proj)
#ax.set_yticks([47.25, 47.5, 47.75, 48], crs=ccrs.Miller())

#ax.add_image(terrain, 9)
#ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'])
#ax.add_feature(rivers_lakes, facecolor=(0.89,0.92,0.94,1), edgecolor='black')
#ax.add_feature(coast, facecolor=(0.89,0.92,0.93,0.8), edgecolor=(0,0,0,0.6))

"""
def plot_sample_grib(gribfiles, title_text="A sample plot"):

    #fig, axs = plt.subplots(2, int(len(gribfiles)/2))
    #
    #for x in range(0, len(gribfiles)):
    #    grb = pygrib.open(gribfiles[x])[1]
    #    ax = axs[x] if len(axs.shape) == 1 else axs
    #    if len(axs.shape) >= 2: ax = axs[int(x/axs.shape[0])][x%axs.shape[0]]

    for f in gribfiles:
        fig, ax = plt.subplots(1, 1)
        grb = pygrib.open(f)[1]

        data = grb.values
        lat, lon = grb.latlons()
        m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='l', ax=ax)
        x, y = m(lon,lat)

        # Paint map with parameter values under projected coordinates.
        #cs = ax.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.jet)
        cs = ax.contourf(x, y, data, cmap=plt.cm.jet)

        # map filigree
        m.drawcoastlines()
        m.fillcontinents()
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
        m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])

        # Plot legend, title.
        fig.colorbar(cs,orientation='vertical', ax=ax)
        #ax.set_title(title_text)
        plt.title(title_text)

        # Show plot.
        fig.tight_layout()
        plt.show()

"""
"""

def plot_coverage(lat, lon):
    fig = plt.figure()
    #m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='c')
    m=Basemap(projection='mill',lat_ts=10,
            llcrnrlon=-180, urcrnrlon=180,
            llcrnrlat=-90,urcrnrlat=90, 
            resolution='c')
    x, y = m(lon,lat)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])
    plt.scatter(x, y, 1, marker='.', color='xkcd:ocean blue', zorder=10)
    fig.tight_layout()
    plt.show()
"""
