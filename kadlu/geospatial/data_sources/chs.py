""" 
    Module within the kadlu package for handling Non-Navigational NONNA-100 
    bathymetric data from the The Canadian Hydrographic Service (CHS). 

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import os
import numpy as np
from kadlu.geospatial.geospatial import crop, read_geotiff
import json
import requests
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg


def fetch(south, north, west, east):
    """ Returns a list of filepaths for downloaded content """

    localfiles = verify_local_files(south, north, west, east) 
    if localfiles is not False: return localfiles

    # api call: get raster IDs within bounding box
    source = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/"
    spatialRel = "esriSpatialRelIntersects"
    spatialReference = "4326"  # WGS-84 spec
    #spatialReference = "3857"  # used for the web map service
    geometry = json.dumps({"xmin":west, "ymin":south, "xmax":east, "ymax":north})
    url1 = f"{source}ImageServer/query?geometry={geometry}&returnIdsOnly=true&geometryType=esriGeometryEnvelope&spatialRel={spatialRel}&f=json&outFields=*&inSR={spatialReference}"
    req1 = requests.get(url1)
    assert(req1.status_code == 200)

    # break up the retrieved objectIds into chunks of size 20
    filepaths = []
    rasterIds = json.loads(req1.text)['objectIds']
    assert(len(rasterIds) > 0)
    for chunk in range(0, int(len(rasterIds) / 20) + 1):
        # api call: get resource location of rasters for each raster ID 
        rasterIdsCSV = ','.join([f"{x}" for x in rasterIds[chunk * 20:(chunk+1) * 20]])
        url2 = f"{source}ImageServer/download?geometry={geometry}&geometryType=esriGeometryPolygon&format=TIFF&f=json&rasterIds={rasterIdsCSV}"
        req2 = requests.get(url2)
        assert(req2.status_code == 200)
        jsondata = json.loads(req2.text)
        assert("error" not in jsondata.keys())

        # api call: for each resource location, download the rasters for each tiff file
        for img in jsondata['rasterFiles']:
            fname = img['id'].split('\\')[-1]
            if 'CA2_' not in fname: continue  # more research required to find out why the Ov_i files exist
            fpath = f"{storage_cfg()}{fname}"
            filepaths.append(fpath)
            if os.path.isfile(fpath):
                print(f"File {fname} exists, skipping retrieval...")
                continue

            print(f"Downloading {fname} from Canadian Hydrographic Service NONNA-100...")
            assert(len(img['rasterIds']) == 1)  # we make an assumption that each file has only one associated raster
            url3 = f"https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/ImageServer/file?id={img['id'][0:]}&rasterId={img['rasterIds'][0]}"
            tiff = requests.get(url3)
            assert(tiff.status_code == 200)
            with open(fpath, "wb") as f: f.write(tiff.content)
        
    return filepaths


def load(south, north, west, east):
    """ 
        Load Non-Navigational NONNA-100 bathymetric data from the Canadian Hydrographic 
        Service (CHS) within specified geographical region.

        Returns 1-dimensional numpy arrays for:
            bathymetry values
            lats
            lons
    """
    files = fetch(south, north, west, east)
    bathy, lats, lons = list(), list(), list()        

    for f in files:
        # load data from a single file
        z,y,x = load_from_file(f) 

        # crop the region of interest
        indices, y, x = crop(y, x, south, north, west, east)
        z = z[indices]

        # collect data arrays
        lats.append(y)
        lons.append(x)
        bathy.append(z)

    # concatenate
    bathy = np.ma.concatenate(bathy)
    lats = np.concatenate(lats)
    lons = np.concatenate(lons)

    return (bathy,lats,lons)


def load_from_file(path):
    """ Load bathymetric data from a GeoTIFF file provided by the Canadian Hydrographic 
        Service (CHS) as part of the Non-Navigational NONNA-100 bathymetric data series.

        Args: 
            path: str
                Path to the GeoTIFF file.

        Returns:
            z: 1d numpy array
                Bathymetry values
            y: 1d numpy array
                Latitude values
            x: 1d numpy array
                Longitude values
    """

    # read data from geotiff file
    #z = read(path)
    z = np.flip(read_geotiff(path=path), axis=0)

    # create lat-lon arrays
    lat, lon = latlon(path)

    # make a grid
    x, y = np.meshgrid(lon, lat)

    # select non-masked entries
    x = x[~z.mask]
    y = y[~z.mask]
    z = z[~z.mask]

    return z,y,x


def latlon(path, num_lat=1001, num_lon=1001):
    """ Create latitude and longitude arrays for a CHS bathymetry data file.

        The number of latitude and longitude values is 1001 by default.

        The latitude binning is 0.001 degrees.

        The longitude binning depends on the latitude, as follows

            * south of 68 deg N: 0.001 deg
            * latitudes from 68 to 80 deg N: 0.002 deg
            * 80 deg N and north: 0.004 deg

        Args: 
            path: str
                File name
            num_lat: int
                Number of latitude grid points
            num_lon: int
                Number of longitude grid points

        Returns:
            lats: numpy array
                Uniformly spaced latitude values 
            lons: numpy array
                Uniformly spaced longitude values 
    """
    # parse SW corner of the map
    south, west = parse_sw_corner(path)

    # latitude step size in degrees
    dlat = 0.001

    # longitude step size in degrees
    if south < 68:
        dlon = 0.001
    elif south >=68 and south < 80:
        dlon = 0.002
    elif south >= 80:
        dlon = 0.004

    lats = np.arange(num_lat, dtype=np.float)
    lats *= dlat
    lats += south

    lons = np.arange(num_lon, dtype=np.float)
    lons *= dlon
    lons += west

    return lats, lons


def parse_sw_corner(path):
    fname = os.path.basename(path)
    south = int(fname[4:8]) / 100
    west = -int(fname[9:14]) / 100
    assert west >= -180 and west <= 180, 'Invalid parsed longitude value'
    assert south >= -90 and south <= 90, 'Invalid parsed latitude value'
    return south, west


def filename(south, west):
    fname = "CA2_{0:04d}N{1:05d}W.tif".format(int(south * 100), -int(west * 100))
    return fname


def verify_local_files(south, north, west, east):
    """ 
        Checks to see if local files exist before querying 
        Returns filenames if true, otherwise returns false.
    """
    fnames = []
    for x in range(int(west), int(east) + 1):
        for y in range(int(south), int(north) + 1):
            f = f"{storage_cfg()}{filename(y, x)}"
            if not os.path.isfile(f): return False
            fnames.append(f)
    print("Files exist, skipping retrieval...")
    return fnames 


class Chs():
    def fetch_bathymetry(self, south=44.4, north=44.7, west=-64.4, east=-63.8):
        return fetch(south, north, west, east)
    def load_bathymetry(self, south=44.4, north=44.7, west=-64.4, east=-63.8):
        return load(south, north, west, east)
    def __str__(self):
        info = "Non-Navigational 100m (NONNA-100) bathymetry dataset from Canadian Hydrographic Datastore"
        args = "(south=-90, north=90, west=-180, east=180)"
        return fetch_util.str_def(self, info, args)

"""
print(Chs())
bathy, lat, lon = Chs().load_bathymetry()
"""

