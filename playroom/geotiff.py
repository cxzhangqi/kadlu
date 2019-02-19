from osgeo import gdal
import numpy as np

path = '/home/oliskir/Meridian/CHS/CHS_NONNA_100/CA2_4300N06000W.tif'

# open dataset
ds = gdal.Open(path)

print(ds.RasterCount)
print(ds.GetMetadata())

band = ds.GetRasterBand(1)

arr = ds.ReadAsArray()
print(arr.shape)


# Get nodata value from the GDAL band object
nodata = band.GetNoDataValue()

#Create a masked array for making calculations without nodata values
arr = np.ma.masked_equal(arr, nodata)
type(arr)

print("[ MAX ] = ", band.GetMaximum())
print("[ SCALE ] = ", band.GetScale())
print("[ UNIT TYPE ] = ", band.GetUnitType())

print(arr.min())

stats = band.GetStatistics( True, True )

print("[ STATS ] =  Minimum={0:.3f}, Maximum={1:.3f}, Mean={2:.3f}, StdDev={3:.3f}".format( \
                stats[0], stats[1], stats[2], stats[3] ))



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm

M = arr.shape[0]
N = arr.shape[1]

X = np.arange(M)
Y = np.arange(N)
X, Y = np.meshgrid(X, Y)

arr = np.flip(arr, axis=0)


fig, ax = plt.subplots()

p = ax.pcolor(X, Y, arr[:M,:N])
cb = fig.colorbar(p)

plt.show()


from rbf.interpolate import RBFInterpolant
import rbf.basis

np.random.seed(1)

basis = rbf.basis.phs2
order = 1

#arr = arr[560:820,660:920]

xx = np.linspace(0, 1, arr.shape[0])
yy = np.linspace(0, 1, arr.shape[1])
xv, yv = np.meshgrid(xx, yy)
xx_obs = np.empty(shape=(arr.shape[0], arr.shape[1], 2))
xx_obs[:,:,0] = xv
xx_obs[:,:,1] = yv
xx_obs = np.reshape(xx_obs, newshape=(arr.shape[0]*arr.shape[1], 2))
uu_obs = np.reshape(arr, newshape=(arr.shape[0]*arr.shape[1]))

x_obs = list()
u_obs = list()
for i in range(xx_obs.shape[0]):
    x = xx_obs[i]
    u = uu_obs[i]
    if u < 0:
        x_obs.append(x)
        u_obs.append(u)

x_obs = np.array(x_obs)
u_obs = np.array(u_obs)

print(x_obs.shape, u_obs.shape)



if False:

    I = RBFInterpolant(x_obs,u_obs, sigma=0.1, basis=basis, order=order)
    vals = np.linspace(0,1,200)
    x_itp = np.reshape(np.meshgrid(vals,vals),(2,200*200)).T # interp points
    u_itp = I(x_itp) # evaluate the interpolant
    # plot the results
    plt.tripcolor(x_itp[:,0],x_itp[:,1],u_itp,vmin=np.min(u_obs),vmax=np.max(u_obs),cmap='viridis')
    plt.scatter(x_obs[:,0],x_obs[:,1],s=100,c=u_obs,vmin=np.min(u_obs),vmax=np.max(u_obs),
                cmap='viridis',edgecolor='k')
    plt.xlim((0.05,0.95))
    plt.ylim((0.05,0.95))
    plt.colorbar()
    plt.tight_layout()
    plt.show()



from scipy.interpolate import griddata



# close dataset
ds = None
band = None


#https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.griddata.html

grid_x, grid_y = np.mgrid[0:1:200j, 0:1:200j]

grid_z0 = griddata(x_obs, u_obs, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(x_obs, u_obs, (grid_x, grid_y), method='linear')
grid_z2 = griddata(x_obs, u_obs, (grid_x, grid_y), method='cubic')

plt.subplot(221)
plt.plot(x_obs[:,0], x_obs[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()