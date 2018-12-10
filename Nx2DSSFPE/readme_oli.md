### environmental data

The bathymetry data is stored in the file `env_database/BathyData_Mariana_500kmx500km.mat`
in the form of three 153 x 145 matrices:

  * latgrat: Latitude at each point of the grid (in degrees)
  * longrat: Longitude at each point of the grid (in degrees)
  * mat: Water depth at each point of the grid (in negative meters)

The rows correspond to fixed latitudes, while the columns correspond 
to fixed longitutedes.

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
