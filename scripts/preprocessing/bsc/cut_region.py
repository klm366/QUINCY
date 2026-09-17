
# Python script to select a region in a netcdf file

# Author: Iria Ayan-Miguez (BSC), 2026

import xarray as xr
import argparse

def read_args():
    parser = argparse.ArgumentParser(description="Script to select the same region than an input netcdf file applied on another netcdf file ")
    parser.add_argument("-r", "--region", type=str, help="Input reference region netcdf file", required=True)
    parser.add_argument("-i", "--input_file", type=str, help="Input file to select region", required=True)
    parser.add_argument("-o", "--output_file", type=str, help="Output file with selected region", required=True)
    
    return parser.parse_args()


if __name__ == "__main__":

  args=read_args()
  region = args.region
  ifile = args.input_file
  ofile = args.output_file

  ds_region = xr.open_dataset(region)
  ds_cut = xr.open_dataset(ifile)

  lon_ds_region_min = ds_region.lon.min().item()
  lon_ds_cut_min = ds_cut.lon.min().item()

  # ds_region [0,360] --> ds_cut [−180,180]
  if lon_ds_region_min >= 0 and lon_ds_cut_min < 0:
    ds_region = ds_region.assign_coords(
        lon=((ds_region.lon + 180) % 360) - 180
    ).sortby("lon")

  # ds_region [−180,180] --> ds_cut [0,360]
  elif lon_ds_region_min < 0 and lon_ds_cut_min >= 0:
    ds_region = ds_region.assign_coords(
        lon=(ds_region.lon % 360)
    ).sortby("lon")

  dx=abs(ds_region.lon.values[1]-ds_region.lon.values[0])
  dy=abs(ds_region.lat.values[1]-ds_region.lat.values[0])

  lon_min = ds_region.lon.min().item()-2*dx
  lon_max = ds_region.lon.max().item()+2*dx
  lat_min = ds_region.lat.min().item()-2*dy
  lat_max = ds_region.lat.max().item()+2*dy
  
  if ds_cut.lat.values[0] > ds_cut.lat.values[-1]:
    lat_reg = slice(lat_max, lat_min)
  else:
    lat_reg = slice(lat_min, lat_max)
  
  subset = ds_cut.sel(
    lon=slice(lon_min, lon_max),
    lat=lat_reg
  )

  subset.to_netcdf(ofile)


