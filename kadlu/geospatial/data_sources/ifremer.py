import os
import sys
import logging
import requests
from ftplib import FTP

import kadlu

def fetch_netcdf_FTP(ftpurl, ftpdir, localdir):
    ftp = FTP(ftpurl)
    ftp.login()
    ftp.cwd(ftpdir)
    with kadlu.Capturing() as output: ftp.retrlines('NLST')
    for ncfile in [f for f in output if kadlu.ext(f, ('.nc',))]:
        logging.info(f'fetching {ncfile}')
        with open(f'{localdir}{ncfile}', 'wb') as fp:
            ftp.retrbinary(f'RETR {ncfile}', fp.write)


class Ifremer():
    def fetch_ifremer_netcdf_hs2013(self): 
        localdir=f'{kadlu.storage_cfg()}testfiles/'
        if not os.path.isdir(localdir): os.mkdir(localdir)
        fetch_netcdf_FTP(ftpurl='ftp.ifremer.fr', ftpdir='/ifremer/ww3/HINDCAST/ATNW/2013_CFSR/hs/', localdir=f'{kadlu.storage_cfg()}testfiles/')


