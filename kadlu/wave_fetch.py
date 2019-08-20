"""
    Wave data fetch module for the kadlu package

    Authors:
        Casey Hilliard      r.casey.hilliard@gmail.com
        Matthew Smith       matthew.smith@dal.ca

    Organization:
        MERIDIAN
        Dalhousie Insitute for Big Data Analytics
     
    Team:
        Acoustic Data Analytics
        Dalhousie University

    Project:
        Kadlu - Tools for underwater soundscape modeling
"""

import numpy as np
from datetime import datetime, timedelta
import os
from os import path
from os.path import dirname
import pygrib
import urllib.request
import configparser
import warnings
import cdsapi
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

waveSources = {
        'ERA5' : {
            'swh' : 'significant_height_of_combined_wind_waves_and_swell',
            'mwd' : 'mean_wave_direction',
            'mwp' : 'mean_wave_period'
        },
        'WWIII' : {
            'swh' : 'hs',
            'mwd' : 'dp',
            'mwp' : 'tp'
        },
        'RDWPS' : {
            'swh' : 'HTSGW',
            'mwd' : 'WVDIR',
            'mwp' : 'PKPER'
        }
    }

regions = {
        'WWIII' : {
            # This may be replaced by a more modular mapping system later
            # depending on how the API works for [region]_[interval] formatting
            'global'    : 'glo_30m',
            'arctic'    : 'ao_30m',
            'pacific'   : 'ep_10m',
            'atlantic'  : {'10m' : 'at_10m', '4m' : 'at_4m'},
            'US west'   : {'10m' : 'wc_10m', '4m' : 'wc_4m'},
            'alaska'    : {'10m' : 'ak_10m', '4m' : 'ak_4m'}
            },
        'RDWPS' : [
            'gulf-st-lawrence',
            'superior',
            'huron-michigan',
            'erie',
            'ontario'
            ]
        }


def init_default_storage_dir(self, msg):
    self.storage_location = (os.path.abspath(dirname(dirname(__file__))) + "/storage/")
    if not os.path.isdir(self.storage_location):
        os.mkdir(self.storage_location)
    warnings.warn("%s storage location will be set to %s" % (msg, self.storage_location))


def instantiate_config(self):
    cfg = configparser.ConfigParser()
    cfg.read(os.path.join(dirname(dirname(__file__)), "config.ini"))
    try:
        self.storage_location = cfg["storage"]["StorageLocation"]
    except KeyError:  # missing config.ini file
        init_default_storage_dir(self, "missing kadlu/config.ini.")

    if self.storage_location is '':  # null value in config.ini
        init_default_storage_dir(self, "null value in kadlu/config.ini.")

    if not os.path.isdir(self.storage_location):  # verify the location exists
        init_default_storage_dir(self, "storage location doesn't exist.")


def validate_wavesource(filepath, sourceDict):
    gribdata = pygrib.fromstring(open(filepath, "rb").read())
    try:
        assert(gribdata['shortName'] in sourceDict.keys() or gribdata['shortName'] in sourceDict.values())
    except AssertionError as err:
        print("Specified source file is not a wave source")
        raise


class WaveFetch():
    def __init__(self):
        instantiate_config(self)

    def fetchERA5(self, wavevar=waveSources['ERA5']['swh'], time=datetime.now()):
        fetchfile = f"{self.storage_location}ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"

        if path.isfile(fetchfile):
            print("File exists, skipping retrieval...")
            validate_wavesource(fetchfile, waveSources['ERA5'])
        else:
            print("Downloading file from Copernicus Climate Data Store...")
            c = cdsapi.Client()
            c.retrieve('reanalysis-era5-single-levels', {
                'product_type'  : 'reanalysis',
                'format'        : 'grib',
                'variable'      : wavevar,
                'year'          : time.strftime("%Y"),
                'month'         : time.strftime("%m"),
                'day'           : time.strftime("%d"),
                'time'          : time.strftime("%H:%M")
            }, fetchfile)

        
    def fetchWWIII(self, wavevar=waveSources['WWIII']['swh'], time=datetime.now(), region=regions['WWIII']['global']):
        fetchname = f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"
        fetchfile = f"{self.storage_location}{fetchname}"

        if path.isfile(fetchfile):
            print("File exists, skipping retrieval...")
            validate_wavesource(fetchfile, waveSources['WWIII'])
        else:
            print("Downloading file from NOAA WaveWatch III...")
            fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fetchname}"
            urllib.request.urlretrieve(fetchurl, fetchfile)


    def fetchRDWPS(self, wavevar=waveSources['RDWPS']['swh'], time=datetime.now(), region=regions['RDWPS'][0]):
        level_ref = 'SFC_0'
        if wavevar is 'WVDIR': level_ref = 'TGL_0'
        hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
        predictionhour = '000'  # any better than 0-hour?
        fetchname = f"CMC_rdwps_{region}_{wavevar}_{level_ref}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"
        fetchfile= f"{self.storage_location}{fetchname}"

        if path.isfile(fetchfile):
            print("File exists, skipping retrieval...")
            validate_wavesource(fetchfile, waveSources['RDWPS'])
        else:
            print("Downloading file from the Regional Deterministic Wave Prediction System...")
            fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/{hour}/{fetchname}"
            urllib.request.urlretrieve(fetchurl, fetchfile)

