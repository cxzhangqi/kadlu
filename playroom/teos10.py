import numpy as np
import gsw
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


#https://teos-10.github.io/GSW-Python/gsw_flat.html#gsw.p_from_z

p = gsw.p_from_z(z=1000, lat=0)
print('p at +1000 m:  {0} dbar'.format(p))
p = gsw.p_from_z(z=-1000, lat=0)
print('p at -1000 m:  {0} dbar'.format(p))
