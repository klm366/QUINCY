# Script to disaggregate CGLS LAI to low/high vegetation cover
# E. Dutra March 2021

from __future__ import print_function
from netCDF4 import Dataset,num2date,date2num
import eccodes as ec
import numpy as np
import pandas as pd
import os
import time 
import sys
import land_cover_veg_utils as pylcu
import xarray as xr

def get_args():
  from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

  description = 'Disagregate CGLS LAI to Low/High vegetation '
  parser = ArgumentParser(description=description,formatter_class=ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i',dest='FINPUTD',type=str,default='NONE',
                     help="Input/Output directory")
  parser.add_argument('-t',dest='FTYPE',type=str,default='grb',
                     help="grb:for grib nc:for netcdf")
  parser.add_argument('-d', dest='YYYY',type=str,default='2019',
                      help="Date (year) of the dataset")
  opts = parser.parse_args()
  return opts  

def read_grib_all(ffile,sname=None):
  
  f = open(ffile,'rb')
  xdata=[]
  while 1:
    gid = ec.codes_grib_new_from_file(f)
    if gid is None: break
    ss = ec.codes_get(gid,'shortName')
    if sname is None:
      xdata.append(ec.codes_get_values(gid))
    elif sname == ss:
      xdata.append(ec.codes_get_values(gid))
    ec.codes_release(gid)
    
  f.close()
  print(f"Read {ss} {len(xdata)} fields from {ffile}")
  return np.array(xdata,dtype='f4')

def read_nc_all(ffile, sname):

  nc = Dataset(ffile, 'r')
  if sname not in nc.variables:
    raise ValueError(f"Variable '{sname}' not found in {ffile}")

  var = nc.variables[sname]
  
  if var.ndim == 3: 
    # (time, lat, lon)
    xdata = [var[it, :, :] for it in range(var.shape[0])]
  elif var.ndim == 2:
    # (lat, lon)
    xdata = [var[:,:]]
  
  nc.close()
  print(f"Read {sname} {len(xdata)} fields from {ffile}")
  return np.ma.array(xdata, dtype='f4')

def fscale_lai(mdata):

  lai={'lai_lv':np.copy(mdata['lai_lv']),
       'lai_hv':np.copy(mdata['lai_hv'])}

  for vty in range(len(lcI['ifs_vtype'])):
    vtyS=lcI['ifs_vtype'][vty].lower()
    if vtyS not in ['l','h']: continue
    tvS="tv"+vtyS;
    cvS="cv"+vtyS;
    laiS='lai_'+vtyS+'v'

    gpT=(mdata['lsm']>0.5) & (mdata[tvS]==lcI['ifs_ind'][vty])
    gpTH=(gpT) & (mdata[cvS]>0.5)
    if np.sum(gpTH) < 100: continue
    lmax=np.max(mdata[laiS],axis=0)
    zmed=np.median(lmax[gpTH])
    print(f"{lcI['ifs_ind'][vty]} {lcI['ifs_name'][vty]} {tvS} {cvS} {np.sum(gpT)} {np.sum(gpTH)} {zmed}")
    alpha=np.minimum(5,np.maximum(0.2,zmed/lmax))
    lai[laiS][:,gpT]=( (1-mdata[cvS][gpT])*alpha[gpT]*mdata[laiS][:,gpT]
                         +mdata[cvS][gpT]*mdata[laiS][:,gpT] )
    
  ## rescale to match aggregated LAI 
  ecveg = mdata['cvl']+mdata['cvh']
  Lecveg=ecveg>0
  laiT = np.where(Lecveg, (mdata['cvl']*lai['lai_lv']+mdata['cvh']*lai['lai_hv'])/ecveg,0)
  zcoef =np.minimum(5,np.maximum(0.2,mdata['lai'][:,Lecveg]/laiT[:,Lecveg]))
  lai['lai_lv'][:,Lecveg] = zcoef*lai['lai_lv'][:,Lecveg]
  lai['lai_hv'][:,Lecveg] = zcoef*lai['lai_hv'][:,Lecveg]
  return lai['lai_lv'],lai['lai_hv'] 

t0=time.time()
opts = get_args()
print(opts)

FINPUTD=opts.FINPUTD
FTYPE=opts.FTYPE
YYYY=opts.YYYY

## Load cross walking table 
lcI = pylcu.load_cross_table_csv(f"{FINPUTD}/esacci_cwt_ifs.csv")

## Read current fields
dlc={}

if FTYPE == 'grb':
  for cvar in ['lsm','cvl','cvh','tvl','tvh']:
    fread=f"{FINPUTD}/{cvar}"
    dlc[cvar] = read_grib_all(f"{FINPUTD}/{cvar}",cvar).squeeze()
  dlc['lai']=read_grib_all(f"{FINPUTD}/month_lai",'lai')
  dlc['lai_lv']=read_grib_all(f"{FINPUTD}/month_lail_v0",'lai_lv')
  dlc['lai_hv']=read_grib_all(f"{FINPUTD}/month_laih_v0",'lai_hv')

elif FTYPE == 'nc':
  for cvar in ['cvl','cvh','tvl','tvh']:
    dlc[cvar] = read_nc_all(f"{FINPUTD}/esacci_lc_agg_ifs_REG_LL.nc",cvar).squeeze()
  for cvar in ['lai', 'lai_lv', 'lai_hv']:
    dlc[cvar]=read_nc_all(f"{FINPUTD}/month_lai_v0",cvar)
  dlc['lsm']=read_nc_all(f"{FINPUTD}/lsm",'lsm').squeeze()
else:
  raise ValueError('FTYPE can only be grb or nc')

## disaggregate LAI and scale to total 
dlc['lai_lv'],dlc['lai_hv'] = fscale_lai(dlc)


## Save fields
if FTYPE == 'grb':
  for cvar,fname in (('lai_lv','lail'),('lai_hv','laih')):
    ftemplate = open(f"{FINPUTD}/month_{fname}_v0",'rb')
    fout = open(f"{FINPUTD}/month_{fname}",'wb')
    ik=0
    while 1:
      gid = ec.codes_grib_new_from_file(ftemplate)
      if gid is None: break
    
      gidC = ec.codes_clone(gid)
      ec.codes_set_values(gidC,dlc[cvar][ik,:])
      ec.codes_write(gidC,fout)
    
      ec.codes_release(gidC)
      ec.codes_release(gid)
      ik=ik+1
    ftemplate.close()
    fout.close()
    print(f"Writen {ik} fields of {cvar} to {fout} using template of {ftemplate}")
 
else:
  fin = f"{FINPUTD}/month_lai_v0"
  fout = f"{FINPUTD}/month_lai_{YYYY}.nc"

  ds = xr.open_dataset(fin, decode_times=False)
  ds = ds.swap_dims({'time': 'month'})

  ds["month"] = xr.DataArray(
    np.arange(0,12,dtype="float32"),
    dims="month",
    attrs={
      "standard_name": "month",
      "long_name": "month",
      "units": f"months since {YYYY}-01-01"
    }
  )

  for cvar, mvar in (('lai_lv', 'Mlail'), ('lai_hv', 'Mlaih'), ('lai', 'Mlai')):
    ds[cvar].values = dlc[cvar]
    ds = ds.rename({cvar: mvar})
  
  ds.to_netcdf(fout,unlimited_dims=[])
  ds.close()

print(f"Finished in {time.time()-t0:.2f} seconds")
