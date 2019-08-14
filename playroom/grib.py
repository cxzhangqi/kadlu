import pygrib

grbs = pygrib.open("../storage/ERA5_reanalysis_significant_height_of_combined_wind_waves_and_swell_2018-01-01_00h.grb2")

for g in grbs:
    print(g)

values = grbs.select(name="Significant height of combined wind waves and swell")[0].values  # numpy array

print(values.shape, values.min(), values.max())

lats,lons = grbs[1].latlons()

print(lats.shape, lons.shape)
