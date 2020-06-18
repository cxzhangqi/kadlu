import os
import logging
import zipfile
import requests

from kadlu.geospatial.data_sources.load_from_file import load_netcdf
from kadlu.geospatial.data_sources.data_util        import          \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized

class Gebco():

    def fetch_bathymetry(self, **kwargs):
        """ fetch gebco netcdf bathymetry, and return the filepath of extracted data """

        if not serialized(seed='gebco_bathy.nc'):
            logging.info('downloading and decompressing gebco bathymetry from netcdf (~8GB)... ')

            # download zipped netcdf to storage directory
            url = 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2020/zip/'
            with requests.get(url, stream=True) as payload_netcdf:
                assert payload_netcdf.status_code == 200, 'error fetching file'
                with open(storage_cfg()+'gebco_netcdf.zip', 'wb') as f:
                    for chunk in payload_netcdf.iter_content(chunk_size=8192): 
                        f.write(chunk)

            # unzip it
            with zipfile.ZipFile(storage_cfg()+'gebco_netcdf.zip', 'r') as zipf:
                zipf.extractall(storage_cfg())
            unzipped = zipfile.ZipFile(storage_cfg()+"gebco_netcdf.zip", "r").namelist()
            ncpath = [_ for _ in unzipped if _[-3:] == '.nc'][0]
            os.rename(storage_cfg() + ncpath, storage_cfg() + 'gebco_bathy.nc')
            renamed = [fname if fname[-3:] != ".nc" else "gebco_bathy.nc" for fname in unzipped]

            # store some metadata
            insert_hash(seed='gebco_bathy.nc')

            logging.info(f'extracted {renamed} to {storage_cfg()}')

        return storage_cfg() + 'gebco_bathy.nc'


    def load_bathymetry(self, **kwargs):
        val, lat, lon = load_netcdf(filename=self.fetch_bathymetry(), **kwargs)
        return val * -1, lat, lon


