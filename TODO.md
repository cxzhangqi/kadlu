### TODO:

- [ ] Exception module
- [ ] Creating analogue to bathy_reader for wave data
    - [ ] wave_fetch -- retrieval of data to "read"
        - [ ] Implement storage configuration
            - [ ] Identify means of user specification of storage space
            - [ ] Reference user specified storage space within fetch functions
                - [ ] Fix function headers referencing storage
        - [ ] Add fetch handler for DalCoast data
        - [ ] Modify all fetch modules to simultaneously fetch 3 parameters (or N parameters as more are integrated into modeling)
        - [ ] Tidy comments / re-check for comments requesting fixes / modifications
        - [ ] Migrate test cases in main () out
    - [ ] wave_reader module
        - [ ] Generalize extraction from GRIB to match with needs of model component (to be provided)
- [ ] Generalize fetch method within wave_fetch 
    - [ ] to handle bathymetry
    - [ ] to handle temperature and salinity
