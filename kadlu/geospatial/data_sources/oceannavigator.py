"""
Reference
https://dfo-ocean-navigator.github.io/Ocean-Navigator-Manual/
"""


import requests
import os
import shutil
from urllib.request import urlopen
from urllib.parse import urlencode
from contextlib import closing
try:
   from PIL import Image
except:
   print("If you are on a Windows machine, please install PIL (imaging library) using 'python -m pip install Pillow' in the Anaconda Prompt")
   exit()
import json

def requestFile():
   # json object for the query
   query = {"area":[{"innerrings":[],
"name":"",
"polygons":[[[45.06576154770312,-61.17187499999999],[43.60426186809619,-64.951171875],[41.153842357114456,-64.79736328124999],[42.33418438593941,-60.84228515624999],[43.779026621608295,-61.61132812499999],[44.535674532413196,-60.46875],[45.06576154770312,-61.17187499999999]]]}],
"bathymetry":1,
"colormap":"default",
"contour":{"colormap":"default",
"hatch":0,
"legend":1,
"levels":"auto",
"variable":"none"},
"dataset":"giops_day",
"depth":"11",
"interp":"gaussian",
"neighbours":10,
"projection":"EPSG:3857",
"quantum":"day",
"quiver":{"colormap":"default",
"magnitude":"length",
"variable":"none"},
"radius":25,
"scale":"30,40,auto",
"showarea":1,
"time":2211451200,
"type":"map",
"variable":"vosaline"}
   # Assemble full request
   base_url = "http://navigator.oceansdata.ca/api/v1.0/plot/?"
   url = base_url + urlencode({"query": json.dumps(query)}) + '&save&format=csv&size=10x7&dpi=144'
   print(url)
   # Save file and finish
   data_file = requests.get(url, stream=True)
   dump = data_file.raw
   # change this if you want a different save location
   location = os.getcwd()
   with open("script_output.csv", "wb") as location:
      print('Saving File')
      shutil.copyfileobj(dump, location)
      print('Done')

if __name__ == '__main__':
   requestFile()
