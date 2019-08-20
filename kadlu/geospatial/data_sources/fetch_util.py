import os
from os import path
from os.path import dirname
import configparser
import warnings
import pygrib

def init_default_storage_dir(msg):
    storage_location = (path.abspath(dirname(dirname(dirname(dirname(__file__))))) + "/storage/")
    if not os.path.isdir(storage_location):
        os.mkdir(storage_location)
    warnings.warn("%s storage location will be set to %s" % (msg, storage_location))
    return storage_location


def instantiate_storage_config():
    cfg = configparser.ConfigParser()
    cfg.read(path.join(dirname(dirname(dirname(dirname(__file__)))), "config.ini"))
    try:
        storage_location = cfg["storage"]["StorageLocation"]
    except KeyError:  # missing config.ini file
        return init_default_storage_dir("missing kadlu/config.ini.")

    if storage_location is '':  # null value in config.ini
        return init_default_storage_dir("null value in kadlu/config.ini.")

    if not path.isdir(storage_location):  # verify the location exists
        return init_default_storage_dir("storage location doesn't exist.")

    return storage_location


def validate_wavesource(filepath, sourceDict):
    gribdata = pygrib.fromstring(open(filepath, "rb").read())
    try:
        assert(gribdata['shortName'] in sourceDict.keys() or gribdata['shortName'] in sourceDict.values())
    except AssertionError as err:
        print("Specified source file is not a wave source")
        raise

