import logging
from PIL import Image
from datetime import datetime

import numpy as np

from kadlu.geospatial.data_sources.data_util        import          \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized


conn, db = database_cfg()

raster_table = lambda var: f'raster_{var}'


def process_rasters_2D(var, filepaths):
    """ process and store arbitrary 2D data from raster format """
    for filepath in filepaths:
        # open image and interpret pixels as elevation
        im = Image.open(filepath)
        nan = float(im.tag[42113][0])
        val = np.ndarray((im.size[0], im.size[1]))
        for yi in range(im.size[1]):
            val[:,yi] = np.array(list(map(im.getpixel, zip(
                    [yi for xi in range(im.size[0])], 
                    range(im.size[1])))))
        mask = np.flip(val == nan, axis=0)

        # generate latlon arrays
        #file_south, file_west = parse_sw_corner(filepath)
        # TODO: 
        # get SW corner of data from file
        assert False, 'need to get SW corner'
        dlat = 0.001
        if file_south < 68:
            dlon = 0.001
        elif file_south >=68 and file_south < 80:
            dlon = 0.002
        elif file_south >= 80:
            dlon = 0.004
        file_xmax = im.size[0] * dlon + file_west
        file_ymax = im.size[1] * dlat + file_south
        file_lon = np.linspace(start=file_west,  stop=file_xmax, num=im.size[0])
        file_lat = np.linspace(start=file_south, stop=file_ymax, num=im.size[1])

        # select non-masked entries, remove missing, build grid
        z1 = np.flip(val, axis=0)
        x1, y1 = np.meshgrid(file_lon, file_lat)
        x2, y2, z2 = x1[~mask], y1[~mask], np.abs(z1[~mask])
        source = ['chs' for z in z2]
        grid = list(map(tuple, np.vstack((z2, y2, x2, source)).T))

        # insert into db
        n1 = db.execute(f"SELECT COUNT(*) FROM {raster_table(var)}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {raster_table(var)} VALUES (?,?,?,?)", grid)
        n2 = db.execute(f"SELECT COUNT(*) FROM {raster_table(var)}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()
        logging.info(f"CHS {filepath.split('/')[-1]} bathymetry in region "
              f"{fmt_coords(dict(south=south,west=west,north=north,east=east))}. "
              f"processed and inserted {n2-n1} rows. "
              f"{len(z1[~mask]) - len(grid)} null values removed, "
              f"{len(grid) - (n2-n1)} duplicate rows ignored")
    
