import numpy as np
from gsw.density import sound_speed

# scalars
c = sound_speed(SA=1, CT=1, p=1)
print(c)

# vectors
SA = np.array([1, 2])
CT = np.array([1, 2])
p = np.array([1, 2])
c = sound_speed(SA=SA, CT=CT, p=p)
print(c)

# multi-d
SA = np.array([[1, 2],[3, 4]])
CT = np.array([[1, 2],[3, 4]])
p = np.array([[1, 2],[3, 4]])
c = sound_speed(SA=SA, CT=CT, p=p)
print(c)
