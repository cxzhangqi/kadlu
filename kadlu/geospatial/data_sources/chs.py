""" 
    Module within the kadlu package for handling Non-Navigational NONNA-100 
    bathymetric data from the The Canadian Hydrographic Service (CHS). 

    Authors: Oliver Kirsebom, Casey Hilliard, Matthew Smith
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
"""

import os
import numpy as np
from kadlu.geospatial.geospatial import crop, read_geotiff
import json
import requests
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from osgeo import gdal


def latlon(filepath):
    """ construct lat / lon arrays based on attributes found in the data """
    bathy = gdal.Open(filepath)
    south, west = parse_sw_corner(filepath)

    dlat = 0.001
    if south < 68:
        dlon = 0.001
    elif south >=68 and south < 80:
        dlon = 0.002
    elif south >= 80:
        dlon = 0.004

    lat = np.arange(start=south, stop=(bathy.RasterYSize*dlat+south), step=dlat, dtype=np.float)
    lon = np.arange(start=west, stop=(bathy.RasterXSize*dlon+west), step=dlon, dtype=np.float)

    return lat, lon


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
    print(storage_cfg())
    # establish which maps are needed
    xmin = int(np.floor(west))
    xmax = int(np.floor(east))
    ymin = int(np.floor(south))
    ymax = int(np.floor(north))
    if abs(east-int(east)) < 1e-6:
        xmax -= 1
    if abs(north-int(north)) < 1e-6:
        ymax -= 1
    # check if maps have already been downloaded
    fnames = []
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            f = os.path.join(storage_cfg(), filename(y, x))
            if not os.path.isfile(f): return False
            fnames.append(f)
    print("Files exist, skipping retrieval...")
    return fnames 


def fetch_chs(south, north, west, east):
    """ Returns a list of filepaths for downloaded content """
    # api call: get raster IDs within bounding box
    #source = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/"
    source = "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/"
    # check source exists
    spatialRel = "esriSpatialRelIntersects"
    spatialReference = "4326"  # WGS-84 spec
    geometry = json.dumps({"xmin":west, "ymin":south, "xmax":east, "ymax":north})
    url1 = f"{source}ImageServer/query?geometry={geometry}&returnIdsOnly=true&geometryType=esriGeometryEnvelope&spatialRel={spatialRel}&f=json&outFields=*&inSR={spatialReference}"
    req1 = requests.get(url1)
    assert(req1.status_code == 200)
    assert("error" not in json.loads(req1.text).keys())

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


def load_chs(south, north, west, east, band_id=1):
    # verify the files exist - if not, fetch them
    fnames = verify_local_files(south, north, west, east)
    if fnames is False: fnames = fetch_chs(south, north, west, east)

    bathy_out = np.array([])
    lat_out = np.array([])
    lon_out = np.array([])
    for filepath in fnames:
        # load data from chs file (bathy,lat,lon)
        z,y,x = load_chs_file(filepath, band_id)

        # discard any data points that are outside requested area
        indices = np.argwhere(np.logical_and(np.logical_and(x >= west, x <= east), np.logical_and(y >= south, y <= north)))
        z = z[indices]
        y = y[indices]
        x = x[indices]

        # append
        bathy_out = np.append(bathy_out, z)
        lat_out = np.append(lat_out, y)
        lon_out = np.append(lon_out, x)

    return bathy_out, lat_out, lon_out


def load_chs_file(filepath, band_id=1):
    # open file
    data_set = gdal.Open(filepath)
    band = data_set.GetRasterBand(band_id)
    values = data_set.ReadAsArray()
    # replace missing values with nan
    nodata = band.GetNoDataValue()
    values[values == nodata] = np.nan
    bathy = np.ma.masked_invalid(values)
    # select non-masked entries
    z = np.flip(bathy, axis=0)
    lat, lon = latlon(filepath)
    x, y = np.meshgrid(lon, lat)
    x = x[~z.mask]
    y = y[~z.mask]
    z = z[~z.mask]

    return z,y,x


class Chs():
    def fetch_bathymetry(self, south=44.4, north=44.7, west=-64.4, east=-63.8):
        return fetch_chs(south, north, west, east)
    def load_bathymetry(self, south=44.4, north=44.7, west=-64.4, east=-63.8):
        return load_chs(south, north, west, east, band_id=1)
    def __str__(self):
        info = "Non-Navigational 100m (NONNA-100) bathymetry dataset from Canadian Hydrographic Datastore"
        args = "(south=-90, north=90, west=-180, east=180)"
        return fetch_util.str_def(self, info, args)

"""
# gulf st lawrence test area:
south =  46
north =  52
west  = -70
east  = -56


print(Chs())
bathy, lat, lon = Chs().load_bathymetry()

filepath = f"{storage_cfg()}CA2_5000N06200W.tif" 
lats, lons = latlon(f"{storage_cfg()}CA2_5000N06200W.tif")
bathy, lat, lon = load_chs(50.5, 50.5, -61.5, -61.5)
"""
