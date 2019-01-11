
from pyost.bathy_reader import BathyReader, LatLon
from pyost.bathy_interpolator import BathyInterpolator
import matplotlib.pyplot as plt

path = '../pyost/tests/assets/BathyData_Mariana_500kmx500km.mat'
reader = BathyReader(path=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat', lon_axis=0)
interp = BathyInterpolator(bathy_reader=reader, rebin_xy=4)

fig = interp.plot_ll()

plt.show()

fig = interp.plot_xy()

plt.show()