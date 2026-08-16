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
bulk_dens=Path(inputDateDir)/'BLDFIE.soilgrid.3600.1800.7.nc'
ph=Path(inputDateDir)/'PHIHOX.soilgrid.3600.1800.7.nc'
pwp_vol=Path(inputDateDir)/'WWP.soilgrid.3600.1800.7.nc'
depth_to_bedrock=Path(inputDateDir)/'BDRICM.soilgrid.3600.1800.nc'
fc_vol=Path(inputDateDir)/'AWCtS.soilgrid.3600.1800.7.nc'

# Open NetCDF files
bulk_dens_dataset=Dataset(bulk_dens,'r',format="NETCDF4")
ph_dataset=Dataset(ph,'r',format="NETCDF4")
pwp_vol_dataset=Dataset(pwp_vol,'r',format="NETCDF4")
depth_to_bedrock_dataset=Dataset(depth_to_bedrock,'r',format="NETCDF4")
fc_vol_dataset=Dataset(fc_vol,'r',format="NETCDF4")

# Extract data from open NetCDF files
BLDFIE=np.array(bulk_dens_dataset['BLDFIE'][:])
BLDFIE_lat=np.array(bulk_dens_dataset['latitude'][:])
BLDFIE_lon=np.array(bulk_dens_dataset['longitude'][:])

PHIHOX=np.array(ph_dataset['PHIHOX'][:])
PHIHOX_lat=np.array(ph_dataset['latitude'][:])
PHIHOX_lon=np.array(ph_dataset['longitude'][:])

WWP=np.array(pwp_vol_dataset['WWP'][:])
WWP_lat=np.array(pwp_vol_dataset['latitude'][:])
WWP_lon=np.array(pwp_vol_dataset['longitude'][:])

BDRICM=np.array(depth_to_bedrock_dataset['BDRICM'][:])

AWCtS=np.array(fc_vol_dataset['AWCtS'][:])
AWCtS_lat=np.array(fc_vol_dataset['latitude'][:])
AWCtS_lon=np.array(fc_vol_dataset['longitude'][:])

# layers 1-5 --> below 1m --> to calculate weighted mean
depth_begin = np.array([0.0, 0.05, 0.15, 0.3, 0.6])
depth_end = np.array([0.05, 0.15, 0.3, 0.6, 1.0])

# Calculate depth_sl (difference between depth_end and depth_begin)
depth_sl = depth_end - depth_begin

# Calculate weighted mean for 5 layers (below 1 meter)
weighted_mean_BLDFIE=np.average(BLDFIE[0:5,:,:],weights=depth_sl,axis=0)
weighted_mean_AWCtS=np.average(AWCtS[0:5,:,:],weights=depth_sl,axis=0)
weighted_mean_WWP=np.average(WWP[0:5,:,:],weights=depth_sl,axis=0)

# Calculate AWCtS
bedrockDepth_cm=BDRICM/100
AWC = (weighted_mean_AWCtS-weighted_mean_WWP)/100*bedrockDepth_cm*1000
AWC[AWC == 0] = -9999

##### Saving files into the 'work_dir' or the pre-processing scripts (var 'outDir') #####

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

# Create a NetCDF file for result_matrix soil pH
create_netcdf_file('PHIHOX.soilgrid.3600.1800.tmp', 1800, 3600, 'PHIHOX', weighted_mean_ph,
                   'Soil pH x 10 in H2O', 'index*10', -9999, PHIHOX_lat, PHIHOX_lon)

# Create a NetCDF file for result_matrix AWCtS / fc_vol
create_netcdf_file('AWCtS.soilgrid.3600.1800.tmp', 1800, 3600, 'AWCtS', AWC,
                   'Saturated water content (volumetric fraction) for tS', 'percent', -9999, WWP_lat, WWP_lon)
