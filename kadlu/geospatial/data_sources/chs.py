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
from osgeo import gdal
import warnings

from kadlu.geospatial.data_sources.data_util import \
storage_cfg, database_cfg, chs_table, str_def


conn, db = database_cfg()


def parse_sw_corner(path):
    """ return the southwest corner coordinates for a given bathymetry file """
    fname = os.path.basename(path)
    south = int(fname[4:8]) / 100
    west = -int(fname[9:14]) / 100
    assert west >= -180 and west <= 180, 'Invalid parsed longitude value'
    assert south >= -90 and south <= 90, 'Invalid parsed latitude value'
    return south, west


# matt_s 2019-12
# removed this function since filenames are now pulled directly from source

#def filename(south, west):
#    """ generate a filename for given southwest corner coordinates """
#    fname = "CA2_{0:04d}N{1:05d}W.tif".format(int(south * 100), -int(west * 100))
#    return fname


def fetch_chs(south, north, west, east, band_id=1):
    """ download bathymetric geotiffs, process them, and insert into db

    args:
        south, north: float
            ymin, ymax coordinate boundaries to fetch bathymetry. range: -90, 90
        west, east: float
            xmin, xmax coordinate boundaries to fetch bathymetry. range: -180, 180

    return: 
      nothing

    """
    # api call: get raster IDs within bounding box
    source = "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/"
    spatialRel = "esriSpatialRelIntersects"
    spatialReference = "4326"  # WGS-84 spec
    geometry = json.dumps({"xmin":west, "ymin":south, "xmax":east, "ymax":north})
    url1 = f"{source}ImageServer/query?geometry={geometry}&returnIdsOnly=true&geometryType=esriGeometryEnvelope&spatialRel={spatialRel}&f=json&outFields=*&inSR={spatialReference}"
    req1 = requests.get(url1)
    assert(req1.status_code == 200)
    assert("error" not in json.loads(req1.text).keys())

    # api call: query for resource locations of rasters
    imgs = []
    rasterIds = json.loads(req1.text)['objectIds']
    assert(len(rasterIds) > 0)
    for chunk in range(0, int(len(rasterIds) / 20) + 1):  # max request size is 20 at a time
        rasterIdsCSV = ','.join([f"{x}" for x in rasterIds[chunk * 20:(chunk+1) * 20]])
        url2 = f"{source}ImageServer/download?geometry={geometry}&geometryType=esriGeometryPolygon&format=TIFF&f=json&rasterIds={rasterIdsCSV}"
        req2 = requests.get(url2)
        assert(req2.status_code == 200)
        jsondata = json.loads(req2.text)
        assert("error" not in jsondata.keys())
        imgs += [img for img in jsondata['rasterFiles'] if 'CA2' in img['id']]

    # api call: for each tiff image, download the associated rasters
    filepaths = []
    imgnum = 1
    print()
    for img in imgs:
        fname = img['id'].split('\\')[-1]
        fpath = f"{storage_cfg()}{fname}"
        print(f"{fname}: downloading {imgnum}/{len(imgs)} from CHS NONNA-100...", end="\r")
        #if os.path.isfile(fpath): continue
        filepaths.append(fpath)
        assert(len(img['rasterIds']) == 1)
        url3 = f"{source}ImageServer/file?id={img['id'][0:]}&rasterId={img['rasterIds'][0]}"
        tiff = requests.get(url3)
        assert(tiff.status_code == 200)
        with open(fpath, "wb") as f: f.write(tiff.content)
        imgnum += 1
        
    print()
    assert(len(filepaths) > 0)

    # read downloaded files and process them for DB insertion
    for filepath in filepaths:
        tiff_data = gdal.Open(filepath)
        band = tiff_data.GetRasterBand(band_id)
        values = tiff_data.ReadAsArray()
        bathy = np.ma.masked_invalid(values)

        # generate latlon arrays
        file_south, file_west = parse_sw_corner(filepath)
        dlat = 0.001
        if file_south < 68:
            dlon = 0.001
        elif file_south >=68 and file_south < 80:
            dlon = 0.002
        elif file_south >= 80:
            dlon = 0.004
        file_ymax = tiff_data.RasterYSize*dlat+file_south
        file_xmax = tiff_data.RasterXSize*dlon+file_west
        file_lat = np.linspace(start=file_south, stop=file_ymax, num=tiff_data.RasterYSize)
        file_lon = np.linspace(start=file_west,  stop=file_xmax, num=tiff_data.RasterXSize)

        # select non-masked entries and remove missing
        z1 = np.flip(bathy, axis=0)
        x1, y1 = np.meshgrid(file_lon, file_lat)
        ix = z1[~z1.mask] != band.GetNoDataValue()
        x2 = x1[~z1.mask][ix]
        y2 = y1[~z1.mask][ix]
        z2 = np.array(z1[~z1.mask][ix].data, dtype=int)

        # build coordinate grid and insert into db
        source = ['chs' for z in z2]
        grid = list(map(tuple, np.vstack((z2, y2, x2, source)).T))
        n1 = db.execute(f"SELECT COUNT(*) FROM {chs_table}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {chs_table} VALUES (?,?,?,?)", grid)
        n2 = db.execute(f"SELECT COUNT(*) FROM {chs_table}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()

        print(f"{filepath.split('/')[-1]} processed and inserted {n2-n1} rows.\t"
              f"{len(z1[~z1.mask]) - len(grid)} null values removed, "
              f"{len(grid) - (n2-n1)} duplicate rows ignored")

    return 


def load_chs(south, north, west, east):
    db.execute(' AND '.join([f"SELECT * FROM {chs_table} WHERE lat >= ?",
                                                              "lat <= ?",
                                                              "lon >= ?",
                                                              "lon <= ?"]),
               tuple(map(str, [south, north, west, east])))
    
    slices = np.array(db.fetchall(), dtype=object).T
    assert len(slices) == 4, "no data found for query range"
    if len(slices) != 4:
       warnings.warn("no data found for query range, returning empty arrays")
       bathy, lat, lon = np.array([]), np.array([]), np.array([])
    else:
       bathy, lat, lon, source = slices
    return np.array((bathy, lat, lon))


class Chs():
    """ collection of module functions for fetching and loading """

    def fetch_bathymetry(self, **kwargs):
        return fetch_chs(south=kwargs['south'], north=kwargs['north'], 
              west=kwargs['west'], east=kwargs['east'], band_id=1)

    def load_bathymetry(self, **kwargs):
        return load_chs(south=kwargs['south'], north=kwargs['north'], 
              west=kwargs['west'], east=kwargs['east'])

    def __str__(self):
        info = "Non-Navigational 100m (NONNA-100) bathymetry dataset from Canadian Hydrographic Datastore"
        args = "(south=-90, north=90, west=-180, east=180)"
        return str_def(self, info, args)


