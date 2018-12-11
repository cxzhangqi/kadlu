### Environmental data

*Bathymetry* 

The bathymetry data are stored in the file `env_database/BathyData_Mariana_500kmx500km.mat`
in the form of three 153 x 145 matrices:

  1. **latgrat**: Latitude at each point of the grid (in degrees)
  2. **longrat**: Longitude at each point of the grid (in degrees)
  3. **mat**: Height of seafloor at each point of the grid (in negative meters)

The rows correspond to fixed latitudes, while the columns correspond 
to fixed longitudes. The binning is roughly 0.03 deg, corresponding 
to a separation of roughly 3 km. The are covered is therefore close 
to 500 x 500 km^2.

*Sound speed profile (SSP)*

The sound speed profile data (along with other environmental variables) 
are stored in the file `WaterColumnProfiles_Mariana_Fake.mat` 
in the form of seven 1 x 54834 matrices (vectors):

  1. **pressuret**: ?
  2. **pressure**: Pressure at given depth (in ?)
  3. **c**: Sound speed at given depth (in meters/second)
  4. **svx**: ?
  5. **Salinity**: Salinity at given depth (in ?)
  6. **Temperature**: Temperature at given depth (in degrees C)
  7. **Depth**: Depth in meters below sea surface (in positive meters)


### create_ENV_models.m 

This script pre-processes the bathymetry and SSP data (see above) and converts 
them into a format suitable for the transmission loss calculation.

*Pre-processing of bathymetry data*

 1. Convert the longitude/latitude values (LL) into position coordinates (XY)
 2. Interpolate the bathymetry grid using the Matlab function `ScatteredInterpolant`.
 3. Create a new XY grid
 4. Compute the bathymetry at each point of the new grid using the interpolator F(X,Y).

*Pre-processing of sound-speed data*

 1. Smoothen the sound-speed profile using a running window of 10 m.
 2. Make a 1D interpolation of smoothened sound speed vs. depth using the Matlab function `interp1`.

The script outputs a file called `Mariana_ENV.mat` containing the two fields `WD` (Water Depth) 
and `SSP` (Sound Speed Profile).


### fRunNx2D_Mariana.m 

This is the main program, which sets up the environment and calculates the 
transmission loss using the Parabolic Equation (PE) solver `propNx2DWAPE`. 
The structure of the program is summarized below.
 
 1. Load the fields `WD` and `SSP`
 2. Specify some additional environmental parameters such as bottom sound speed, 
    bottom attenuation, and bottom density.
 3. Specify source depth and frequency.
 4. Specify a number of settings for the PE solver. 
 5. Run the PE solver.
 6. Save the results and make some nice plots.


### subroutines/propNx2DWAPE.m

This is the PE solver.


### calculate_noise_fields.m

This script does some post-processing of the output of the main program.
