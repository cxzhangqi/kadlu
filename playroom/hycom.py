"""
Experimenting with HYCOM - OPeNDAP Dataset Access

Access domain:
    https://tds.hycom.org/thredds/dodsC/
Example data access web API: 
    https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
Example data response page (JSON):
    https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.ascii?depth%5B0:1:39%5D,lat%5B0:1:3250%5D,lon%5B0:1:4499%5D,time%5B0:1:2860%5D,water_temp%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,water_temp_bottom%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,salinity%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,salinity_bottom%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,surf_el%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D
"""

import numpy as np
import urllib
import json

url = r'https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.ascii?depth%5B0:1:39%5D,lat%5B0:1:3250%5D,lon%5B0:1:4499%5D,time%5B0:1:2860%5D,water_temp%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,water_temp_bottom%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,salinity%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,salinity_bottom%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D,surf_el%5B0:1:0%5D%5B0:1:0%5D%5B0:1:0%5D'

jsonstring = urllib.request.urlopen(url)
data = jsonstring.read().decode('utf8')
if jsonstring.msg == 200:  # status ok
    jdata = json.loads(data)

#binarydata = open("playroom/2015.dods", "rb").read()
