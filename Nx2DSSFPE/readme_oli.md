### Environmental data

*Bathymetry* 

The bathymetry data are stored in the file `env_database/BathyData_Mariana_500kmx500km.mat`
in the form of three 153 x 145 matrices:

  1. latgrat: Latitude at each point of the grid (in degrees)
  2. longrat: Longitude at each point of the grid (in degrees)
  3. mat: Water depth at each point of the grid (in negative meters)

The rows correspond to fixed latitudes, while the columns correspond 
to fixed longitutedes.

*Sound speed profile (SSP)*

The data needed to compute the sound speed are stored in the file `WaterColumnProfiles_Mariana_Fake.mat` 
in the form of seven 1 x 54834 matrices (vectors):

  1. pressuret
  2. pressure
  3. = c
  4. = svx
  5. = Salinity
  6. = Temperature
  7. = Depth


### create_ENV_models.m 

This script takes bathymetry data and sound-speed profile (SSP) data from 
the `env_database` folder as input and outputs `Mariana_ENV.mat`, which serves 
as input for the main program.


### fRunNx2D_Mariana.m 

This is the main program, which sets up the environment and calculates the 
transmission loss using the PE solver.


### subroutines/propNx2DWAPE.m

This is the PE solver script.


### calculate_noise_fields.m

This script does some post-processing of the output of the main program.
