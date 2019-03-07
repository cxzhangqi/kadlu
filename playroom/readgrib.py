
# Source for GRIB files of wave, wind, precip, etc. 
# spanning St. Lawrence:
# https://weather.gc.ca/grib/grib2_Gulf-St-Lawrence_e.html


import numpy as np
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime

### ECMWF ERA5 Format:
# Source: see ecmwf_era5_grib_fetch.py, cdsapi API calls
def fetchERA5():

    # Sample target file
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/03_Sample_ECMWF_ERA5/era5_reanalysis_sig_wave_swell_2018_Jan_00_hourly.grib'

    # Fetch grib structure from target.
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

    # Attempt fetch of matching data.
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

def fetchWWIII():
    # NOAA WAVEWatch III Format:
    # Source: URL form https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/YYYY/MM/gribs/multi_1.glo_30m.hs.WWWWMM.grb2 
    # e.g. URL https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/2018/11/gribs/multi_1.glo_30m.hs.201811.grb2


    # Sample target file
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201811.grb2'
    #grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201712.grb2'
    #grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/01_Sample_Wavewatch_III/multi_1.glo_30m.hs.201810.grb2'

    # Fetch grib structure from target.
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

def fetchRDWPS():
    
    # CMC RDWPS (Gulf of St. Lawrence, regional forecast)
    # Source (daily only): URL: http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/HH/CMC_rdwps_gulf-st-lawrence_PARAM_SFC_0_latlon0.05x0.05_YYYYMMDD00_PTTT.grib2
    # HH Forecast product hour 00/03/06/09; PARAM - forecast parameter; YYYYMMDD (Current date); TTT - Target forecast out in 3hr intervals from 000 (?nowcast) to 48 hrs out
    grib = '/home/hilliard/MARIN/08_BigData/50_MERIDIAN/05_kadlu/02_Sample_CMC_StL/CMC_rdwps_gulf-st-lawrence_HTSGW_SFC_0_latlon0.05x0.05_2019022100_P000.grib2'

    # Fetch grib structure from target.
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

    # Fetch parameter values, coordinates.
    data=grb.values
    lat,lon = grb.latlons()

    # Fetch, project basemap.
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

    # Test ERA5
#    (grbsample, samp_title_text) = fetchERA5()
#    plotSampleGrib(grbsample, samp_title_text)
    
    # Test WWIII
#    (grbsample, samp_title_text) = fetchWWIII()
#    plotSampleGrib(grbsample, samp_title_text)
    
    # Test RDWPS
    (grbsample, samp_title_text) = fetchRDWPS()
    plotSampleGrib(grbsample, samp_title_text)
    

if __name__ == '__main__':
    main()

