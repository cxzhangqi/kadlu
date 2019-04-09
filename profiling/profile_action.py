import cProfile
from kadlu.transmission_loss_calculator import TransmissionLossCalculator

calculator = TransmissionLossCalculator(bathymetry=None, sound_speed=None, step_size=None, range=50e3, angular_bin=10, vertical_bin=10)

cProfile.run('calculator.run(frequency=10, source_depth=9900)', 'calc.prof')

# Usage:
#  python profile_action.py
#  snakeviz calc.prof