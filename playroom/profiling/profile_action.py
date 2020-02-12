import sys
import cProfile
import numpy as np
from kadlu.transmission_loss.transmission_loss_calculator import TransmissionLossCalculator
from kadlu.geospatial.data_provider import DataProvider

n = int(sys.argv[1])

print('\n(Remember to recompile and reinstall Kadlu before profiling)')
print('Profiling case #{0} ...'.format(n))

if n == 0:
    calculator = TransmissionLossCalculator(env_data=None, sound_speed=None, flat_seafloor_depth=10000, step_size=None, range=50e3, angular_bin_size=5, vertical_bin_size=10, steps_btw_bathy_updates=1)
    cProfile.run('calculator.run(frequency=10, source_depth=9900)', 'calc.prof')

elif n == 1:
    folder = "../kadlu/tests/assets/tif"
    provider = DataProvider(storage_location=folder, bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    seafloor_depth = -provider.bathy(x=0, y=0)
    max_depth = -np.min(provider.bathy_data[0]) 
    calculator = TransmissionLossCalculator(env_data=provider, sound_speed=None,\
        step_size=None, range=30e3, angular_bin_size=5, vertical_bin_size=10,\
        max_depth=1.2*max_depth, steps_btw_sound_speed_updates=1,\
        verbose=False, progress_bar=True)
    cProfile.run('calculator.run(frequency=10, source_depth=0.9*seafloor_depth)', 'calc.prof')


print()

# Usage:
#  python profile_action.py
#  snakeviz calc.prof
