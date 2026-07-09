# Script to aggregate CGLS LAI or FCOVER to a coarser grid 
# E. Dutra  March 2021

from __future__ import print_function
from netCDF4 import Dataset,num2date,date2num
import numpy as np
import time 
import datetime as dt
import eccodes as ec
import os 

import land_cover_veg_utils as pylcu
import land_cover_veg_cython_utils as cext


def gen_nc(fname,gout,cvarN,fclim=True):
  
  print("Generating:",fname)
  os.system('rm -f '+fname)
  
  ncfile = Dataset(fname,mode='w',format='NETCDF4')
  ncfile.history = 'Created ' + time.ctime(time.time())
  
  ## dimension
  if not fclim:
    ncfile.createDimension('time',None)
  else:
    ncfile.createDimension('time',12)
    
  if gout['nlon'] == 0 : 
    ncfile.createDimension('x',gout['nptot'])
    ncfile.createDimension('y',1)
    chk3=[1,1,gout['nptot']]
    dim3=('time','y','x')
  
    ## Variables
    cvar=ncfile.createVariable('lat','f8',('y','x'))
    cvar.long_name = 'latitude'
    cvar.standard_name = 'latitude'
    cvar.units = 'degrees_north'
    cvar.axis = "Y"
    cvar[:] = gout['lat']
    
    cvar=ncfile.createVariable('lon','f8',('y','x'))
    cvar.long_name = 'longitude'
    cvar.standard_name = 'longitude'
    cvar.units = 'degrees_east'
    cvar.axis = "X"
    cvar[:] = gout['lon']

  else:
    ncfile.createDimension('lon',gout['nlon'])
    ncfile.createDimension('lat',gout['nlat'])
    chk3=[1,gout['nlat'],gout['nlon']]
    dim3=('time','lat','lon')
    
    cvar=ncfile.createVariable('lat','f8',('lat'))
    cvar.long_name = 'latitude'
    cvar.standard_name = 'latitude'
    cvar.units = 'degrees_north'
    cvar._CoordinateAxisType = "Lat"
    cvar[:] = gout['latr']
    
    cvar=ncfile.createVariable('lon','f8',('lon'))
    cvar.long_name = 'longitude'
    cvar.standard_name = 'longitude'
    cvar.units = 'degrees_east'
    cvar._CoordinateAxisType = "Lon" 
    cvar[:] = gout['lonr']
  
  cvar=ncfile.createVariable('time','f8',('time'))
  cvar.long_name = 'Time'
  cvar.calendar = 'standard'
  cvar.units = 'days since 1970-01-01 00:00:00'
  cvar[0] = 0
  if fclim:
    cdates = [ date2num(dt.datetime(1970,im,1),cvar.units) for im in range(1,13) ]
    cvar[:] = cdates
  if cvarN=='LAI':
    cvars=['lai_lv','lai_hv','lai']
  for cv in cvars:
    cvar=ncfile.createVariable(cv,'f4',dim3,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk3)
    cvar.long_name = cvarN
    cvar.units=''
    cvar.coordinates = "lat lon"

  
  
  return ncfile
  
def get_args():
  from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

  description = 'Aggregate CGLS LAI/FCOVER '
  parser = ArgumentParser(description=description,formatter_class=ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i',dest='FINPUT',type=str,default='NONE',
                     help="Input file")
  parser.add_argument('-cvar',dest='CVAR',type=str,default='LAI',
                     help="Variable (LAI or FCOVER)")
  parser.add_argument('-o',dest='FOUT',type=str,default='NONE',
                     help="output file")
  parser.add_argument('-s',dest='subarea',type=int,default=0,
                     help="0 for global, 1 for domain")
  parser.add_argument('-cci',dest='FCCI',type=str,default='NONE',
                     help="high resolutio land cover matching input")
  parser.add_argument('-g',dest='FGOUT',type=str,default='NONE',
                     help="File with target grid for output")
  parser.add_argument('-d',dest='sdate',type=str,default='NONE',
                     help="date in format YYYYMMDD")
  opts = parser.parse_args()
  return opts  

t0=time.time()
opts = get_args()
print(opts)

CVAR=opts.CVAR
FOUT=opts.FOUT
FINPUT=opts.FINPUT
FCCI=opts.FCCI
FGOUT=opts.FGOUT
subarea=opts.subarea
SDATE=opts.sdate
STEP_READ=0

zmiss=np.double(-999.)
WDATE=dt.datetime.strptime(SDATE,'%Y%m%d')


## Read input / output grids 
gIN = pylcu.get_grid_reg(FINPUT)
gOUT = pylcu.get_grid_reg(FGOUT)
if gOUT['nlon']==0:
  LGG=True
else:
  LGG=False

## Load high resolution land cover 
t = time.time()
nc = Dataset(FCCI,'r')
lsm = nc.variables['lsm'][:]
cvl = nc.variables['cvl'][:]
cvh = nc.variables['cvh'][:]
nc.close()
print("Loading HRES land cover %s %5.1f sec"%(FCCI,time.time()-t))

## Open input file & read data
t = time.time()
nc = Dataset(FINPUT,'r')
data_in = nc.variables[CVAR][STEP_READ,:,:]
data_in = np.ma.filled(data_in,zmiss)
nc.close()
print("Loading LAI data %s in %5.1f sec"%(FINPUT,time.time()-t))

# Allocate output array 
# 3 dimensions mlail,mlaih,mlai
if LGG:
  gIN['lonr']=pylcu.wrapTo360(gIN['lonr'])
  data_out = np.zeros((3,1,gOUT['nptot']),dtype='f4')
else:
  data_out = np.zeros((3,gOUT['nlat'],gOUT['nlon']),dtype='f4')

## Do the aggregation    
t=time.time()
if LGG:
  cext.agg_cgls_red(gIN['latr'],gIN['lonr'],
                    data_in,lsm,cvl,cvh,data_out,
                    gOUT['latr'].astype('f4'),
                    gOUT['nplon'].astype('i4'),
                    gOUT['nplonC'],
                    zmiss)
else: 
  cext.agg_cgls_reg(gIN['latr'],gIN['lonr'],
                    data_in,lsm,cvl,cvh,data_out,
                    gOUT['latr'][0],gOUT['dlat'],
                    gOUT['lonr'][0],gOUT['dlon'],
                    gOUT['nlat'],gOUT['nlon'],
                    subarea,zmiss)
    
data_out = np.ma.masked_where(data_out[:,:,:]==zmiss,data_out[:,:,:])
  
print("Aggregating in %5.1f sec"%(time.time()-t))
  
## Save to output
t=time.time()

if LGG:
  ftemplate = open(FGOUT,'rb')
  gid_template = ec.codes_grib_new_from_file(ftemplate)
  ## ftemplate must contain LSM
  sn = ec.codes_get(gid_template,'shortName')
  if sn != 'lsm':
    raise ValueError('grib template provide is not LSM!')
  ilsm = ec.codes_get_values(gid_template)
  isnotland = ilsm <= 0.5
  idate=int(WDATE.strftime("%Y%m%d"))
  foutgrb = open(FOUT,'wb')
  
  for ik,pparam in enumerate([66,67,500206]):
    xout = data_out[ik,:].filled(0).squeeze()
    xout[isnotland] = 0 
    kk={'param':pparam,"dataDate":idate,}
    if pparam == 500206:  ec.codes_set(gid_template,'centre','cnmc')
    pylcu.save2grb(xout,gid_template,foutgrb,kk)

  
  ftemplate.close()
  foutgrb.close()
  
else:
  print("Creating:",FOUT)
  ncOUT=gen_nc(FOUT,gOUT,CVAR,fclim=False)
  ncOUT.variables['lai_lv'][0,:,:]=data_out[0,:,:]
  ncOUT.variables['lai_hv'][0,:,:]=data_out[1,:,:]
  ncOUT.variables['lai'][0,:,:]=data_out[2,:,:]
  ncOUT.variables['time'][0] = date2num(WDATE,ncOUT.variables['time'].units)
  ncOUT.close()
  
print("Saving  in %5.1f sec"%(time.time()-t))

print("Finished in in %5.1f sec"%(time.time()-t0))
