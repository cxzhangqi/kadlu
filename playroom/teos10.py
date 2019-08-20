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



# complete example
# starting point is practical salinity (SP) and in-situ temperature (t)

# depths
z = np.array([-10, -130, -3200])

# latitudes
lats = np.array([43.0, 44.1, 45.2])

# longitudes
lons = np.array([-60.1, -60.2, -60.3])

# practical salinities
SP = np.array([35.1, 35.0, 34.9])

# in-situ temperatures
t = np.array([4.1, 4.2, 4.3])

# compute sea pressure
p = gsw.p_from_z(z=z, lat=lats)

# compute absolute salinities
SA = gsw.SA_from_SP(SP, p, lons, lats)

# compute conservative temperatures
CT = gsw.CT_from_t(SA, t, p)

# compute sound speed
c = sound_speed(SA=SA, CT=CT, p=p)

print('sea pressures: ',p)
print('absolute salinities: ',SA)
print('conservatime temperatures: ',CT)
print('sound speeds: ', c)
