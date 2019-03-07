
# Source for GRIB files of wave, wind, precip, etc. 
# spanning St. Lawrence:
# https://weather.gc.ca/grib/grib2_Gulf-St-Lawrence_e.html


import numpy as np
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime
from enum import Enum
import urllib.request

def loadERA5():
### ECMWF ERA5 Format:
# Source: see ecmwf_era5_grib_load.py, cdsapi API calls

    # Sample target file
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/03_Sample_ECMWF_ERA5/era5_reanalysis_sig_wave_swell_2018_Jan_00_hourly.grib'

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Identify parameter and date for extraction.
    # Date field: validDate
    # Target date (not incl. time yet)
    date_valid = datetime(2018,1,29)

    # Short ECMWF codes:
    #   swh - Sig wave height
    #   mwd - Mean wave direction (degrees)
    #   mwp - Mean wave period (s)
    # Short ECMWF field: shortNameECMF --OR-- shortName 

    # Attempt load of matching data.
    #grb = grbs.select(validDate=date_valid,name='Significant height of combined wind waves and swell')[0]
    grb = grbs.select(validDate=date_valid,shortNameECMF='swh')[0]

    # Build title text to match.
    title_text = "Example 3: ECMWF ERA5 Wave Height from GRIB\n({}) ".format(date_valid)

    # Spatial components
    # longitudeOfFirstGridPointInDegrees
    # latitutdeOfFirstGridPointInDegrees
    # Ni
    # Nj
    return (grb, title_text)

# Enumerated class defining WWIII Reanalysis coverages.
class WWIIIRegion(Enum):
    glo_30m = "Global 30 min"
    ao_30m = "Arctic Ocean 30 min"
    at_10m = "NW Atlantic 10 min"
    wc_10m = "US West Coast 10 min"
    ep_10m = "East Pacific 10 min"
    ak_10m = "Alaskan 10 min"
    at_4m = "Gulf of Mexico and NW Atlantic 4 min"
    wc_4m = "US West Coast 4 min"
    ak_4m = "Alaskan 4 min"

# Enumerated class defining WWIII Reanalysis parameters to use.
class WWIIIWavevar(Enum):
    hs = "swh"
    dp = "mwd"
    tp = "mwp"
    
def fetchWWIII(fetch_timestamp=datetime.now(), region=WWIIIRegion.glo_30m, wavevar=WWIIIWavevar.hs):

    # Peel off strings from fetch date for component parts of the fetchc url
    fetch_year = fetch_timestamp.strftime("%Y") 
    fetch_month = fetch_timestamp.strftime("%m")
    
    # Build URL
    fetchURLString = 'https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/' + fetch_year + '/' + fetch_month + '/gribs/multi_1.' + region.name + '.' + wavevar.name + '.' + fetch_year + fetch_month + '.grb2'
    fetchURLFile = 'multi_1.' + region.name + '.' + wavevar.name + '.' + fetch_year + fetch_month + '.grb2'

    ### DEBUG ###
    print("URL:{}\nFILE:{}".format(fetchURLString, fetchURLFile))

    # Attempt to retrieve the referneced target.
    urllib.request.urlretrieve(fetchURLString, fetchURLFile)

    ### Handle storage here.

def loadWWIII():
    # NOAA WAVEWatch III Format:
    # Source: URL form https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/YYYY/MM/gribs/multi_1.glo_30m.hs.WWWWMM.grb2 
    # e.g. URL https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/2018/11/gribs/multi_1.glo_30m.hs.201811.grb2


    # Sample target file
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201811.grb2'
    #grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201712.grb2'
    #grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201810.grb2'

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Identify parameter and date for extraction.
    # Date field: validDate
    # Target date (not incl. time yet)
    date_valid = datetime(2018,11,2)

    # Short ECMWF codes are used here too, WWIII tokens appear to exist only on the file wrapper (see line grib = ...)
    #   swh - Sig wave height
    #   mwd - Mean wave direction (degrees)
    #   mwp - Mean wave period (s)
    # Short ECMWF field: shortNameECMF --OR-- shortName 
    grb = grbs.select(validDate=date_valid,shortNameECMF='swh')[0]

    title_text = "Example 2: NWW3 Sig. Wave Height from GRIB\n({}) ".format(date_valid)

    return (grb, title_text)

# Enumerated class defining RDWPS coverages.
class RDWPSRegion(Enum):
    gulf_st_lawrence = "Gulf of St. Lawrence"
    superior = "Lake Superior"
    huron_michigan = "Lake Huron and Michigan"
    erie = "Lake Erie"
    ontario = "Lake Ontario"

# Enumerated class defining RDWPS parameters to use.
class RDWPSWavevar(Enum):
    HTSGW = "swh"
    WVDIR = "mwd" ## Check this, other sources do not refer to only wind driven for direction (CH 20190307).
    PKPER = "mwp"
    
def fetchRDWPS(fetch_timestamp=datetime.now(), region=RDWPSRegion.gulf_st_lawrence, wavevar=RDWPSWavevar.HTSGW):

    # Peel off strings from fetch date for component parts of the fetchc url
    fetch_year = fetch_timestamp.strftime("%Y") 
    fetch_month = fetch_timestamp.strftime("%m")
    fetch_day = fetch_timestamp.strftime("%d")
    
    ## If we can constrain to: 00, 03, 06, 09, we can allow specification
    # -- are any better than 00 for our use?
    fetch_forecast_hour = '00'
    
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
    
    # Build URL
    fetchURLString = 'http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/' + fetch_forecast_hour + '/CMC_rdwps_' + hyph_region_name + '_' + wavevar.name + '_' + level_ref + '_latlon0.05x0.05_' + fetch_year + fetch_month + fetch_day + '00_P' + fetch_prediction_hour + '.grib2'
    fetchURLFile = 'CMC_rdwps_' + hyph_region_name + '_' + wavevar.name + '_' + level_ref + '_latlon0.05x0.05_' + fetch_year + fetch_month + fetch_day + '00_P' + fetch_prediction_hour + '.grib2'

    ### DEBUG ###
    print("URL:{}\nFILE:{}".format(fetchURLString, fetchURLFile))

    # Attempt to retrieve the referneced target.
    urllib.request.urlretrieve(fetchURLString, fetchURLFile)

    ### Handle storage here.


def loadRDWPS():
    
    # CMC RDWPS (Gulf of St. Lawrence, regional forecast)
    # Source (daily only): URL: http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/HH/CMC_rdwps_gulf-st-lawrence_PARAM_SFC_0_latlon0.05x0.05_YYYYMMDD00_PTTT.grib2
    # HH Forecast product hour 00/03/06/09; PARAM - forecast parameter; YYYYMMDD (Current date); TTT - Target forecast out in 3hr intervals from 000 (?nowcast) to 48 hrs out
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/02_Sample_CMC_StL/CMC_rdwps_gulf-st-lawrence_HTSGW_SFC_0_latlon0.05x0.05_2019022100_P000.grib2'

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Identify parameter and date for extraction.
    # Date field: validDate
    # Target date (not incl. time yet)
    date_valid = datetime(2019,2,21)

    #hs	- significant height of combined wind waves and swell, metres (HTSGW)
    #tp	- primary wave mean period, seconds (PKPER)
    #dp	- primary wind wave direction, degrees true (i.e. 0 deg = North; proceeding clockwise) (WVDIR)
    # Short ECMWF codes are used here too, RDWPS tokens appear to exist only on the file wrapper (see line grib = ...)
    #   swh - Sig wave height
    #   mwd - Mean wave direction (degrees)
    #   mwp - Mean wave period (s)
    #grb = grbs.select()[0]
    grb = grbs.select(validDate=date_valid,shortNameECMF='swh')[0]    
    
    title_text = "Example 1: CMC-RDWPS GoSL Sig. Wave Height from GRIB\n({}) ".format(date_valid)

    return (grb, title_text)

def plotSampleGrib(grb,title_text):

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
    m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

    # Plot legend, title.
    plt.colorbar(cs,orientation='vertical')
    plt.title(title_text)

    # Show plot.
    plt.show()

def main():

    # Test WWIII Fetch
#    fetchWWIII(datetime(2017,2,3))
    # Test RDWPS Fetch
#    fetchRDWPS()

    # Test ERA5 Load
#    (grbsample, samp_title_text) = loadERA5()
#    plotSampleGrib(grbsample, samp_title_text)
    
    # Test WWIII Load
#    (grbsample, samp_title_text) = loadWWIII()
#    plotSampleGrib(grbsample, samp_title_text)
    
    # Test RDWPS Load
    (grbsample, samp_title_text) = loadRDWPS()
    plotSampleGrib(grbsample, samp_title_text)
    

if __name__ == '__main__':
    main()

