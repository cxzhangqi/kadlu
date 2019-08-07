#!/usr/bin/env python3
""" Wave data fetch module within the kadlu package

    Authors: Casey Hilliard
    contact: r.casey.hilliard@gmail.com
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from enum import Enum
import urllib.request
import os.path
from os.path import dirname
import configparser
import warnings

# API for fetching from ECMWF data, requires definition of .cdsapirc file with 
# download URL / key pair.
import cdsapi

class WaveSources(Enum):
    """ Enum class for supported wave sources."""
    ECMWF_ERA5 = 0
    NOAA_WWIII = 1
    EC_RDWPS = 2

class ERA5Wavevar(Enum):
    """ Enum class for wave parameter names in ERA5 source
    """

    significant_height_of_combined_wind_waves_and_swell = "swh"
    mean_wave_direction = "mwd"
    mean_wave_period = "mwp"

class WWIIIWavevar(Enum):
    """ Enum class for wave parameter names in NOAA WAVEWATCH III source
    """

    hs = "swh"
    dp = "mwd"
    tp = "mwp"

class RDWPSWavevar(Enum):
    """ Enum class for wave parameter names in EC RDWPS source
    """

    HTSGW = "swh"
## Check WVDIR, other sources do not refer to only wind driven for direction (CH 20190307).
    WVDIR = "mwd" 
    PKPER = "mwp"


class RDWPSRegion(Enum):
    """ Enum class forecast region names in EC RDWPS source
    """

    gulf_st_lawrence = "Gulf of St. Lawrence"
    superior = "Lake Superior"
    huron_michigan = "Lake Huron and Michigan"
    erie = "Lake Erie"
    ontario = "Lake Ontario"

class WWIIIRegion(Enum):
    """ Enum class forecast region names in NOAA WAVEWATCH III source
    """

    glo_30m = "Global 30 min"
    ao_30m = "Arctic Ocean 30 min"
    at_10m = "NW Atlantic 10 min"
    wc_10m = "US West Coast 10 min"
    ep_10m = "East Pacific 10 min"
    ak_10m = "Alaskan 10 min"
    at_4m = "Gulf of Mexico and NW Atlantic 4 min"
    wc_4m = "US West Coast 4 min"
    ak_4m = "Alaskan 4 min"


class WaveFetch():
    """ Class for fetching / locating wave data.
    
        WaveFetch is used to locate GRIB-based model output wave data, generally
        to be handed off to a WaveReader-like module. Data files can either be
        read from local storage, or fetched from one of several web sources.
        WaveReader expects input GRIB files to contain, at minimum, parameters
        for 1) Significant combined wave and swell height, 2) Mean wave 
        direction in degrees and 3) Mean wave period. Parameter naming is 
        expected to adhere to the 

        Attributes: 
            storage_location: str
                Normally, this should be the path where either a) data exist to
                be loaded, or b) where data should be written to disk.
            fetch_date: Timestamp
                Timestamp for which data are sought.
            wave_source: str
                Enum () value corresponding to the source to be utilized. 
                See ().
            
    """
    def __init__(self, storage_location=None, fetch_datetimestamp=datetime.now(), wave_source=WaveSources.ECMWF_ERA5):

### ADD VALIDATION (file -- exists and source -- isa -> WaveSource) (CH 20190418)

        def init_default_storage_dir(self, msg):
            self.storage_location = (os.path.abspath(dirname(dirname(__file__))) + "/storage/")
            if not os.path.isdir(self.storage_location):
                os.mkdir(self.storage_location)
            warnings.warn("%s storage location will be set to %s" % (msg, self.storage_location))

        self.storage_location = storage_location

        if self.storage_location is None:  # by default. read from config.ini 
            cfg = configparser.ConfigParser()
            cfg.read(os.path.join(dirname(dirname(__file__)), "config.ini"))
            try:
                self.storage_location = cfg["storage"]["StorageLocation"]
            except KeyError:  # missing config file
                init_default_storage_dir(self, "missing kadlu/config.ini.")
        elif self.storage_location is '':  # null value in config.ini
            init_default_storage_dir(self, "null value in kadlu/config.ini.")

        if not path.isdir(self.storage_location):  # verify the location exists
            init_default_storage_dir(self, "storage location doesn't exist.")

        # TODO:
        # ensure the source is a WaveSource
  
        # Date and source.
        self.fetch_datetimestamp = fetch_datetimestamp
        self.wave_source = wave_source

    
    # Linux/UNIX config: https://cds.climate.copernicus.eu/api-how-to
    # Windows config: https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows
    def fetchERA5(self, wavevar=ERA5Wavevar.significant_height_of_combined_wind_waves_and_swell,fetch_timestamp=None):
        """ Locate and download grib file of ERA5 modeled wave parameter data to 
            storage location, given the variable of interest and timestamp for 
            the model run. 

            Args: 
                wavevar: ERA5Wavevar
                    Variable to be fetched, ERA5 nomenclature.
                fetch_timestamp: datetime
                    Date of file to be fetched.

            Returns:
                None.
        """

        # Use module default timestamp if not overridden.
        if(fetch_timestamp is None):
            fetch_timestamp = self.fetch_datetimestamp

        # Instantiate the cdsapi client.
        c = cdsapi.Client()

        # Peel off strings from fetch date for component parts of cdsapi request.
        fetch_year = fetch_timestamp.strftime("%Y") 
        fetch_month = fetch_timestamp.strftime("%m")
        fetch_day = fetch_timestamp.strftime("%d")
        fetch_hour = fetch_timestamp.strftime("%H:%M")

        # Build a target filename under which the data will be stored.
        self.fetch_filename = self.storage_location + 'ERA5_reanalysis_{}_{}.grb2'.format(wavevar.name, fetch_timestamp.strftime("%Y-%m-%d_%Hh"))

        # Establish a request to obtain the targeted wavevar at the timestamp indicated.
        # Check if a file under the target name already exists, abort fetch if so.
        if os.path.isfile(self.fetch_filename):
            print("File exists, skipping retrieval.")

        # Attempt to retrieve the referenced target, if necessary.
        else:

            c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type':'reanalysis',
                'format':'grib',
                'variable':wavevar.name,
                'year':fetch_year,
                'month':fetch_month,
                'day':fetch_day,
                'time':fetch_hour
            },
            self.fetch_filename)
        
    def loadERA5(self, target_date = None, wavevar=ERA5Wavevar.significant_height_of_combined_wind_waves_and_swell, grib=None):
        """ Loads a single time interval slice from a specified grib file of 
            ERA5 modeled wave parameter data. Returns the grib structure and 
            descriptive text indicating the extraction extent. 

            Args: 
                target_date: datetime
                    Date of internal data product to be fetched.
                wavevar: ERA5Wavevar
                    Variable to be loaded.
                grib: string
                    Path and filename of grib file from which extraction is to 
                    be loaded.

            Returns:
                grb: grib structure object
                    A single variable slice of grib data to plot (data, lat, 
                    lon).
                title_text:
                    Title text snippet to be used in plotting, indicating the 
                    variable and time slice.
        """

        # If no gribfile argument is provided, default to the fetched file.
        if grib is None:
            grib = self.fetch_filename

        # load grib structure from target.
        grbs=pygrib.open(grib)

        # Identify parameter and date for extraction.
        # Date field: validDate
        # Target date (not incl. time yet)
        if (target_date is None):
            target_date = self.fetch_datetimestamp.replace(minute=0, hour=0, second=0, microsecond=0)
        else:
### Should any filtering on time be added here to enforce valid time intervals?
            target_date = target_date

        # Fetch the indicated slice from the overall Grib file.
        grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]

        if(wavevar == ERA5Wavevar.significant_height_of_combined_wind_waves_and_swell):
            title_text = "ERA5 Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
        elif(wavevar == ERA5Wavevar.mean_wave_period):
            title_text = "ERA5 Mean Wave Period from GRIB\n({}) ".format(target_date)
        elif(wavevar == ERA5Wavevar.mean_wave_direction):
            title_text = "ERA5 Mean Wave Direction from GRIB\n({}) ".format(target_date)
        else:
            title_text = "ERA5\nUnknown variable from GRIB\n({}) ".format(target_date)

        return (grb, title_text)

        
    def fetchWWIII(self, region=WWIIIRegion.glo_30m, wavevar=WWIIIWavevar.hs, fetch_timestamp=None):
        """ Locate and download grib file of WWIII modelled wave parameter data to 
            storage location, given the WWIII model region, variable of 
            interest and timestamp for the model run. 

            Args: 
                region: WWIIIRegion
                    Region among those modeled by NOAA as WWIII.
                wavevar: WWIIIWavevar
                    Variable to be fetched, WWIII nomenclature.
                fetch_timestamp: datetime
                    Date of file to be fetched.

            Returns:
                None.
        """

        # Use module default timestamp if not overridden.
        if(fetch_timestamp is None):
            fetch_timestamp = self.fetch_datetimestamp

        # Peel off strings from fetch date for component parts of the fetchc url
        fetch_year = fetch_timestamp.strftime("%Y") 
        fetch_month = fetch_timestamp.strftime("%m")
        
        # Build URL, target filename.
        fetchURLString = 'https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/' + fetch_year + '/' + fetch_month + '/gribs/multi_1.' + region.name + '.' + wavevar.name + '.' + fetch_year + fetch_month + '.grb2'
        self.fetch_filename = self.storage_location + 'multi_1.' + region.name + '.' + wavevar.name + '.' + fetch_year + fetch_month + '.grb2'

        # Check if a file under the target name already exists, abort fetch if so.
        if os.path.isfile(self.fetch_filename):
            print("File exists, skipping retrieval.")

        # Attempt to retrieve the referenced target, if necessary.
        else:
            
            urllib.request.urlretrieve(fetchURLString, self.fetch_filename)

    def loadWWIII(self, target_date = None, wavevar=WWIIIWavevar.hs, grib=None):
        """ Loads a single time interval slice from a specified grib file of 
            WWIII modeled wave parameter data. Returns the grib structure and 
            descriptive text indicating the extraction extent. 

            Args: 
                target_date: datetime
                    Date of internal data product to be fetched.
                wavevar: WWIIIWavevar
                    Variable to be loaded.
                grib: string
                    Path and filename of grib file from which extraction is to 
                    be loaded.

            Returns:
                grb: grib structure object
                    A single variable slice of grib data to plot (data, lat, 
                    lon).
                title_text:
                    Title text snippet to be used in plotting, indicating the 
                    variable and time slice.
        """

        # If no gribfile argument is provided, default to the fetched file.
        if grib is None:
            grib = self.fetch_filename

        # load grib structure from target.
        grbs=pygrib.open(grib)

        # Identify parameter and date for extraction.
        # Date field: validDate
        # Target date (not incl. time yet)
# DEBUG        date_valid = datetime(2018,11,2)
        if (target_date is None):
            target_date = self.fetch_datetimestamp.replace(minute=0, hour=0, second=0, microsecond=0)
        else:
### Should any filtering on time be added here to enforce valid time intervals?
            target_date = target_date

        # Fetch the indicated slice from the overall Grib file.
        grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]

        if(wavevar == WWIIIWavevar.hs):
            title_text = "WW III Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
        elif(wavevar == WWIIIWavevar.tp):
            title_text = "WW III Mean Wave Period from GRIB\n({}) ".format(target_date)
        elif(wavevar == WWIIIWavevar.dp):
            title_text = "WW III Mean Wave Direction from GRIB\n({}) ".format(target_date)
        else:
            title_text = "WW III\nUnknown variable from GRIB\n({}) ".format(target_date)

        return (grb, title_text)

    def fetchRDWPS(self, region=RDWPSRegion.gulf_st_lawrence, wavevar=RDWPSWavevar.HTSGW, fetch_timestamp=None):
        """ Locate and download grib file of EC/CMC modelled wave parameter data 
            to storage location, given the RDWPS model region, variable
            of interest and timestamp for the model run. 

            Args: 
                region: RDWPSRegion
                    Region among those modeled by EC / CMC as RDWPS.
                wavevar: RDWPSWavevar
                    Variable to be fetched, EC nomenclature.
                fetch_timestamp: datetime
                    Date of file to be fetched.

            Returns:
                None.
        """

        # Use module default timestamp if not overridden.
        if(fetch_timestamp is None):
            fetch_timestamp = self.fetch_datetimestamp

        # Peel off strings from fetch date for component parts of the fetchc url
        fetch_year = fetch_timestamp.strftime("%Y") 
        fetch_month = fetch_timestamp.strftime("%m")
        fetch_day = fetch_timestamp.strftime("%d")
        
        ## Truncate the fetch hour to nearest earlier 6-hour interval (00, 06, 12, 18)
        # -- is there a better option? 
        fetch_forecast_hour = "{:02d}".format(((fetch_timestamp.hour % 24) // 6) * 6)
        
        ## If we can constrain to: 000 -> 048, in intervals of 3, we can allow specification
        # -- are any better than 0-hour for our use?
        fetch_prediction_hour = '000'

        # Un-sanitize Enum name, re-add hyphen for url building.
        hyph_region_name = region.name.replace('_','-')

        # Build appropriate level reference based on parameter.
        if (wavevar is RDWPSWavevar.HTSGW or wavevar is RDWPSWavevar.PKPER):
            level_ref = 'SFC_0'
        else: # For WVDIR
            level_ref = 'TGL_0'
        
        # Build URL targeting CMC RDWPS (In test: Gulf of St. Lawrence, regional forecast)
        # Source (daily only): URL: http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/HH/CMC_rdwps_gulf-st-lawrence_PARAM_SFC_0_latlon0.05x0.05_YYYYMMDD00_PTTT.grib2
        # HH Forecast product hour 00/03/06/09; PARAM - forecast parameter; YYYYMMDD (Current date); TTT - Target forecast out in 3hr intervals from 000 (?nowcast) to 48 hrs out
        fetchURLString = 'http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/' + fetch_forecast_hour + '/CMC_rdwps_' + hyph_region_name + '_' + wavevar.name + '_' + level_ref + '_latlon0.05x0.05_' + fetch_year + fetch_month + fetch_day + fetch_forecast_hour + '_P' + fetch_prediction_hour + '.grib2'
        self.fetch_filename = self.storage_location + 'CMC_rdwps_' + hyph_region_name + '_' + wavevar.name + '_' + level_ref + '_latlon0.05x0.05_' + fetch_year + fetch_month + fetch_day + fetch_forecast_hour + '_P' + fetch_prediction_hour + '.grib2'

        # Check if a file under the target name already exists, abort fetch if so.
        if os.path.isfile(self.fetch_filename):
            print("File exists, skipping retrieval.")

        # Attempt to retrieve the referenced target, if necessary.
        else:

            urllib.request.urlretrieve(fetchURLString, self.fetch_filename)

    def loadRDWPS(self, target_date = None, wavevar=RDWPSWavevar.HTSGW, grib=None):
        """ Loads a single time interval slice from a specified grib file of 
            modeleted RDWPS wave parameter data. Returns the grib structure and 
            descriptive text indicating the extraction extent. 

            Args: 
                target_date: datetime
                    Date of internal data product to be fetched.
                wavevar: RDWPSWavevar
                    Variable to be loaded.
                grib: string
                    Path and filename of grib file from which extraction is to be 
                    loaded.

            Returns:
                grb: grib structure object
                    A single variable slice of grib data to plot (data, lat, lon).                    
                title_text:
                    Title text snippet to be used in plotting, indicating the 
                    variable and time slice.
        """

        # target_date: date for extraction.
        # grib: target path/file for extraction

        # If no gribfile argument is provided, default to the fetched file.
        if grib is None:
            grib = self.fetch_filename

        # If no date argument is provided, default to the module timestamp, truncated to nearest earlier 6-hour interval.
        if target_date is None:

            rep_hour = ((self.fetch_datetimestamp.hour % 24) // 6) * 6
            target_date = self.fetch_datetimestamp.replace(minute=0, hour=rep_hour, second=0, microsecond=0)
        else:
### Should similar filtering on time be added here to enforce 6 hour intervals?
            target_date = target_date

        # load grib structure from target.
        grbs=pygrib.open(grib)

        # Fetch slice of data corresponding to selected target date and wave variable.
            #swh	- significant height of combined wind waves and swell, metres (HTSGW)
            #mwp	- primary wave mean period, seconds (PKPER)
            #mwd	- primary wind wave direction, degrees true (i.e. 0 deg = North; proceeding clockwise) (WVDIR)
        grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]
        
        if(wavevar == RDWPSWavevar.HTSGW):
            title_text = "CMC-RDWPS Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
        elif(wavevar == RDWPSWavevar.PKPER):
            title_text = "CMC-RDWPS Mean Wave Period from GRIB\n({}) ".format(target_date)
        elif(wavevar == RDWPSWavevar.WVDIR):
            title_text = "CMC-RDWPS Mean Wave Direction from GRIB\n({}) ".format(target_date)
        else:
            title_text = "CMC-RDWPS\nUnknown variable from GRIB\n({}) ".format(target_date)

        return (grb, title_text)

    def plotSampleGrib(self,grb,title_text):
        """ Plots a single time interval slice from a specified grib file to
        test functionality. Takes the gribfile slice and title text as arguments,
        generates a matplotlib plot as output. 

            Args: 
                grb: grib structure object
                    A single variable slice of grib data to plot (data, lat, lon).
                title_text: string
                    String to appear as title text for the plot.

            Returns:
                None.
        """

        # Instantiate plot
        plt.figure()

        # load parameter values, coordinates.
        data=grb.values
        lat,lon = grb.latlons()

        # load, project basemap.
        m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='l')
          
        # Project data coordinates
        x, y = m(lon,lat)

        # Paint map with parameter values under projected coordinates.
        cs = m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.jet)

        # Draw and fill basemap boundary, map filigree.
        m.drawcoastlines()
        m.fillcontinents()
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
        m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])

        # Plot legend, title.
        plt.colorbar(cs,orientation='vertical')
        plt.title(title_text)

        # Show plot.
        plt.show()

def main():

    # Test ERA5 Fetch
    wfTestERA5 = WaveFetch('/home/hilliard/storage/kadlu_storage/', datetime(2018, 1, 1, 0, 0, 0, 0), WaveSources.ECMWF_ERA5)
    wfTestERA5.fetchERA5()
    # Test ERA5 Load
    (grbsample, samp_title_text) = wfTestERA5.loadERA5()
    wfTestERA5.plotSampleGrib(grbsample, samp_title_text)
    '''    
    # Test WWIII end to end
#    wfTestWWIII = WaveFetch('/tmp/', datetime(2017, 2, 3, 0, 0, 0, 0), WaveSources.ECMWF_ERA5)
    wfTestWWIII = WaveFetch('/home/hilliard/storage/kadlu_storage/', datetime(2017, 2, 3, 0, 0, 0, 0), WaveSources.NOAA_WWIII)
    wfTestWWIII.fetchWWIII()
    (grbsampleWWIII, samp_title_text_WWIII) = wfTestWWIII.loadWWIII()
    wfTestWWIII.plotSampleGrib(grbsampleWWIII, samp_title_text_WWIII)
    
    # Test RDWPS end to end
#    wfTestRDWPS = WaveFetch('/tmp/', datetime.now()  - timedelta(hours=3), WaveSources.EC_RDWPS)
    wfTestRDWPS = WaveFetch('/home/hilliard/storage/kadlu_storage/', datetime.now()  - timedelta(hours=3), WaveSources.EC_RDWPS)
    wfTestRDWPS.fetchRDWPS()
    (grbsampleRDWPS, samp_title_textRDWPS) = wfTestRDWPS.loadRDWPS()
    wfTestRDWPS.plotSampleGrib(grbsampleRDWPS, samp_title_textRDWPS)
    '''
if __name__ == '__main__':
    main()

