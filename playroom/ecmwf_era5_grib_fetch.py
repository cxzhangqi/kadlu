#!/usr/bin/env python
import cdsapi

# Sample fetch
#c = cdsapi.Client()
#c.retrieve("reanalysis-era5-pressure-levels",
#{
#"variable": "temperature",
#"pressure_level": "1000",
#"product_type": "reanalysis",
#"year": "2008",
#"month": "01",
#"day": "01",
#"time": "12:00",
#"format": "grib"
#},
#"download.grib")


### Retrieval for ERA reanalysis: significant wave+swell height, 2018, all days, all times
'''c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'grib',
        'variable':'significant_height_of_combined_wind_waves_and_swell',
        'year':'2018',
        'month':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ]
    },
    'era5_reanalysis_sig_wave_swell_2018_hourly.grib')
'''

### Retrieval for ERA reanalysis: significant wave+swell height, 2018, all days, 0h time
c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'grib',
        'variable':'significant_height_of_combined_wind_waves_and_swell',
        'year':'2018',
        'month':[
            '01'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00'
        ]
    },
    'era5_reanalysis_sig_wave_swell_2018_Jan_00_hourly.grib')
