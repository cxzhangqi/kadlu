import os
import sys
import logging
import requests
from ftplib import FTP

import kadlu


class Ifremer():

    def fetch_hs():
        ftp = FTP('ftp.ifremer.fr')
        ftp.login()
        ftp.cwd('/ifremer/ww3/HINDCAST/ATNW/2013_CFSR/hs/')
        with kadlu.Capturing() as output: ftp.retrlines('NLST')
        for ncfile in [f for f in output if kadlu.ext(f, ('.nc',))]:
            logging.info(f'fetching {ncfile}')
            with open(f'{kadlu.storage_cfg()}testfiles/{ncfile}', 'wb') as fp:
                ftp.retrbinary(f'RETR {ncfile}', fp.write)

