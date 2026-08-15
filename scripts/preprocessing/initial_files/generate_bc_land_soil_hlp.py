# Import necessary libraries
from netCDF4 import Dataset
import warnings
import math
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib as mpl
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
from pathlib import Path
import sys
import os
warnings.filterwarnings("ignore")

# get two arguments, absolute path of input and putput directory
inputDateDir=sys.argv[1]
outDir=sys.argv[2]

# Input files
bulk_dens=inputDateDir/"BLDFIE.soilgrid.3600.1800.7.nc"

# Open NetCDF files
bulk_dens_dataset=Dataset(bulk_dens,'r',format="NETCDF4")

# Extract data from open NetCDF files
BLDFIE=np.array(bulk_dens_dataset['BLDFIE'][:])
BLDFIE_lat=np.array(bulk_dens_dataset['latitude'][:])
BLDFIE_lon=np.array(bulk_dens_dataset['longitude'][:])

# layers 1-5 --> below 1m --> to calculate weighted mean
depth_begin = np.array([0.0, 0.05, 0.15, 0.3, 0.6])
depth_end = np.array([0.05, 0.15, 0.3, 0.6, 1.0])

# Calculate depth_sl (difference between depth_end and depth_begin)
depth_sl = depth_end - depth_begin
depth_sl

# Calculate weighted mean for 5 layers (below 1 meter)
weighted_mean_BLDFIE=np.average(BLDFIE[0:5,:,:],weights=depth_sl,axis=0)


##### Saving files #####

# creating function creating NetCDF file
def create_netcdf_file(file_path, lat_dim_size, lon_dim_size, variable_name, data, long_name, units, fill_value, lat_values, lon_values):
    os.chdir(outDir)
    ncfile = Dataset(file_path, 'w', format="NETCDF4")

    lat_dim = ncfile.createDimension('lat', lat_dim_size)
    lat = ncfile.createVariable('lat', np.double, ('lat',))
    lat[:]=lat_values
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat.standard_name='latitude'
    lat.axis='Y'
    lat.missing_value = -9999
    lat.scale_factor = 1
    lat.add_offset = 0

    lon_dim = ncfile.createDimension('lon', lon_dim_size)
    lon = ncfile.createVariable('lon', np.double, ('lon',))
    lon[:]=lon_values
    lon.standard_name = "longitude"
    lon.long_name = "longitude"
    lon.units = "degrees_east"
    lon.axis = "X"
    lon.missing_value = -9999
    lon.scale_factor = 1
    lon.add_offset = 0

    variable = ncfile.createVariable(variable_name, np.float32, ('lat', 'lon'), fill_value=fill_value)
    ncfile[variable_name][:] = data
    variable.long_name = long_name
    variable.units = units
    variable.missing_value = fill_value

    ncfile.close()


# Create a NetCDF file for result_matrix bulk density
create_netcdf_file('BLDFIE_3600_1800.tmp', 1800, 3600, 'BLDFIE', weighted_mean_BLDFIE,
                   'Bulk density (fine earth) in kg / cubic-meter', 'kg / cubic-m', -9999, BLDFIE_lat, BLDFIE_lon)

