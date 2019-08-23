import cProfile
from kadlu.transmission_loss.transmission_loss_calculator import TransmissionLossCalculator

calculator = TransmissionLossCalculator(env_data=None, sound_speed=None, flat_seafloor_depth=10000, step_size=None, range=50e3, angular_bin_size=1, vertical_bin_size=10)

cProfile.run('calculator.run(frequency=10, source_depth=9900, vertical_slice=False)', 'calc.prof')

# Usage:
#  python profile_action.py
#  snakeviz calc.prof