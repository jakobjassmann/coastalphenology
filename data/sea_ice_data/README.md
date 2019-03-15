# Sea-Ice Data Readme

These scripts prepare NOAA/NSIDC Passive Microwave Sea Ice Concentration Data available through the NSIDC via https://nsidc.org/data/G02202 on a per site basis, extracting mean daily sea-ice extent within a approx 500 km x 500 km bounding box around each site.

'''
/data/sea_ice_data/alexfiord/sea_ice_extraction_alexfiord.R
/data/sea_ice_data/barrow/sea_ice_extraction_barrow.R
/data/sea_ice_data/qhi/sea_ice_extraction_qhi.R
/data/sea_ice_data/zackenberg/sea_ice_extraction_zackenberg.R

'''

The NetCDF files for the time-period covered by the phenology datasets (1990-2016) will have to be downloaded locally, and the references to the NetCDF file locations need to be adjusted accordingly within each script. 
