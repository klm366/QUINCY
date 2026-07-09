# Script to aggregate ESA CCI land cover to a coarser grid 
# E. Dutra March 2021

from __future__ import print_function
from netCDF4 import Dataset,num2date,date2num
import numpy as np
import os
import time 
import sys

import land_cover_veg_utils as pylcu
import land_cover_veg_cython_utils as cext



def gen_nc_ifs(fname,grid,lcI):
  
  print("Generating:",fname)
  os.system('rm -f '+fname)
  
  ncfile = Dataset(fname,mode='w',format='NETCDF4')
  ncfile.history = 'Created ' + time.ctime(time.time())
  
  ## dimension
  if grid['nlon'] == 0 : 
    ncfile.createDimension('x',gout['nptot'])
    ncfile.createDimension('y',1)
    chk3=[1,1,np.minimum(gout['nptot'],chunkS**2)]
    dim3=('vtype','y','x')
    chk2=[1,np.minimum(gout['nptot'],chunkS**2)]
    dim2=('y','x')
    
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
    chk3=[1,1,gout['nlon']]
    dim3=('vtype','lat','lon')
    chk2=[1,gout['nlon']]
    dim2=('lat','lon')
    
    cvar=ncfile.createVariable('lat','f8',('lat'))
    cvar.long_name = 'latitude'
    cvar.standard_name = 'latitude'
    cvar.units = 'degrees_north'
    cvar._CoordinateAxisType = "Lat"
    cvar[:] = grid['latr']
    
    cvar=ncfile.createVariable('lon','f8',('lon'))
    cvar.long_name = 'longitude'
    cvar.standard_name = 'longitude'
    cvar.units = 'degrees_east'
    cvar._CoordinateAxisType = "Lon" 
    cvar[:] = grid['lonr']
  
  cvar=ncfile.createVariable('cvl','f4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk2)
  cvar.long_name = 'low vegetation cover'
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  cvar=ncfile.createVariable('cvh','f4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk2)
  cvar.long_name = 'high vegetation cover'
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  cvar=ncfile.createVariable('tvh','f4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk2)
  cvar.long_name = 'dominant high vegetation type'
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  cvar=ncfile.createVariable('tvl','f4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk2)
  cvar.long_name = 'dominant low vegetation type'
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  cvar=ncfile.createVariable('lsm','f4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk2)
  cvar.long_name = 'land sea (water) mask '
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  return ncfile 

def gen_nc_agg(fname,gout,lcI):
  
  print("Generating:",fname)
  os.system('rm -f '+fname)
  
  ncfile = Dataset(fname,mode='w',format='NETCDF4')
  ncfile.history = 'Created ' + time.ctime(time.time())
  
  ## dimension
  ncfile.createDimension('class',lcI['ncci'])
  if gout['nlon'] == 0 : 
    ncfile.createDimension('x',gout['nptot'])
    ncfile.createDimension('y',1)
    chk=[1,1,np.minimum(gout['nptot'],chunkS**2)]
    dim3=('class','y','x')
    dim2=('y','x')
    
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
    chk=[1,np.minimum(gout['nlat'],chunkS),np.minimum(gout['nlon'],chunkS)]
    dim3=('class','lat','lon')
    dim2=('lat','lon')
    
    ## Variables
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
    cvar._CoordinateAxisType = "Lat"
    cvar[:] = gout['lonr']
  
  cvar=ncfile.createVariable('class','i4',('class'))
  cvar.long_name = 'cci_class_number'
  cvar.standard_name = 'cci_class_number'
  cvar.units = '-'
  cvar.axis = "Z"
  cvar[:] = np.arange(0,lcI['ncci'])
  
  
  cvar=ncfile.createVariable('class_name',str,('class'))
  cvar.long_name = 'cci_class_name'
  cvar.standard_name = 'cci_class_name'
  cvar.units = '-'
  for k in range(lcI['ncci']):
    cvar[k] = "%02i %s"%(k,str(lcI['cci_name'][k]))
  
  cvar=ncfile.createVariable('classfr','f4',dim3,
                             fill_value=zmiss,zlib=True,complevel=1,
                             chunksizes=chk)
  cvar.long_name = 'cci_class_fraction'
  cvar.standard_name = 'cci_class_fraction'
  cvar.units = '-'
  cvar.coordinates = "lat lon"
  
  cvar=ncfile.createVariable('totcount','i4',dim2,
                             fill_value=zmiss,zlib=True,complevel=1)
  cvar.long_name = 'total_pixels_cci'
  cvar.standard_name = 'total_pixels_cci'
  cvar.units = '-'
  cvar.coordinates = "lat lon"

  
  return ncfile   

def get_args():
  from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

  description = 'Aggregate ESA-CCI Land Cover '
  parser = ArgumentParser(description=description,formatter_class=ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i',dest='FCCI',type=str,default='NONE',
                     help="Input CCI file")
  parser.add_argument('-c',dest='FLCCT2PFS',type=str,default='NONE',
                     help="csv cross walking table")
  parser.add_argument('-g',dest='FGOUT',type=str,default='NONE',
                     help="File with grid for output")
  parser.add_argument('-s',dest='subarea',type=int,default=0,
                     help="0 for global, 1 for domain")
  parser.add_argument('-o',dest='FOUTP',type=str,default='NONE',
                     help="Output directory path")
  parser.add_argument('-t',dest='FOUTT',type=str,default='NONE',
                     help="Output tag")
  opts = parser.parse_args()
  return opts  

zmiss=-9999.
chunkS=1024
fmin=0.001  # minimum land franction for calculations 

opts = get_args()
print(opts)
FCCI = opts.FCCI 
FLCCT2PFS=opts.FLCCT2PFS
FGOUT=opts.FGOUT
subarea=opts.subarea
FOUTP=opts.FOUTP
FOUTT=opts.FOUTT

FOUTAGG="%s/esacci_lc_agg_%s.nc"%(FOUTP,FOUTT)
FOUTIFS="%s/esacci_lc_agg_ifs_%s.nc"%(FOUTP,FOUTT)

LAGG=True
LIFS=True 
tstart0=time.time()



if LAGG:
  tstart=time.time()

  ## load cci info and cross walking table 
  lcI = pylcu.load_cross_table_csv(FLCCT2PFS)
  ## create look-up table for the aggregation 
  lcI['cci_indI']=np.zeros(np.max(lcI['cci_ind'])+1,dtype='uint16')+np.max(lcI['cci_ind'])+1
  for ik in range(lcI['ncci']):
    lcI['cci_indI'][lcI['cci_ind'][ik]]=ik
  
  ## get output grid description
  gout = pylcu.get_grid_reg(FGOUT)
  if gout['nlon']==0:
    LGG=True
  else:
    LGG=False
    
  ## get CCI grid description
  gcci = pylcu.get_grid_reg(FCCI)

  ###==================================================
  ### Process in slots of the same size of the netcdf chucks to optimize read speed 
  #### create output array to call cython function
  
  if LGG:
    gcci['lonr']=pylcu.wrapTo360(gcci['lonr'])
    cci_count = np.zeros((lcI['ncci'],1,gout['nptot']),dtype='uint16')
  else:
    cci_count = np.zeros((lcI['ncci'],gout['nlat'],gout['nlon']),dtype='uint16')

  nc_cci = Dataset(FCCI,'r')
  chunks= nc_cci.variables['lccs_class'].chunking()
  lat_chk = chunks[1]*4
  lon_chk = chunks[2]*4
  
  ik=0
  t0=time.time()
  nslots=len(np.arange(0,gcci['nlat'],lat_chk))*len(np.arange(0,gcci['nlon'],lon_chk))
  print('LAGG chunks:',lat_chk,lon_chk,nslots)
  for ilat in np.arange(0,gcci['nlat'],lat_chk):
    for ilon in np.arange(0,gcci['nlon'],lon_chk):
      t = time.time()
      xcci_raw = nc_cci.variables['lccs_class'][0,ilat:ilat+lat_chk,ilon:ilon+lon_chk]
      if LGG : 
        cext.agg_cci_red(gcci['latr'][ilat:ilat+lat_chk],
                     gcci['lonr'][ilon:ilon+lon_chk],
                     xcci_raw,
                     cci_count,
                     lcI['cci_indI'],
                     gout['latr'].astype('f4'),gout['nplon'].astype('i4'),gout['nplonC'])     
      else: 
      
        cext.agg_cci_reg(gcci['latr'][ilat:ilat+lat_chk],
                     gcci['lonr'][ilon:ilon+lon_chk],
                     xcci_raw,
                     cci_count,
                     lcI['cci_indI'],
                     gout['latr'][0],gout['dlat'],
                     gout['lonr'][0],gout['dlon'],
                     gout['nlat'],gout['nlon'],
                     subarea)
      ik=ik+1
      
      sys.stdout.write("\rProcessing :%i of %i chunks [%3i %%] %f sec"%(ik,nslots,100.*ik/nslots,time.time()-t));sys.stdout.flush()  
  nc_cci.close()
  
  ## Normalize
  cci_Ncount = np.sum(cci_count,axis=0,dtype='f4')
  
  print(' ')
  print("Processed CCI aggregation in %6.2f seconds"%(time.time()-t0))

  if not LIFS or gout['nlon'] == 1440:
    ## Only write this intermediate step if LIFS is not active or for regular lat/lon
    ## Generate the output file
    t0 = time.time()
    nc = gen_nc_agg(FOUTAGG,gout,lcI)
      
    ## Save to the output 
    nc.variables['totcount'][:]=cci_Ncount
    for k in range(lcI['ncci']):
      sys.stdout.write("\rSaving field %i"%k);sys.stdout.flush()    
      nc.variables['classfr'][k,:,:] = np.ma.fix_invalid(cci_count[k,:,:]/cci_Ncount)
    nc.close()
    print(' ')
    print("Saving classfr data in %6.2f seconds"%(time.time()-t0))
  
  
  print("Finished LAGG in %6.2f seconds"%(time.time()-tstart))


if LIFS:
  ## Aggregate cci:
  tstart=time.time()

  if not LAGG:
    ## load cci info and cross walking table 
    lcI = pylcu.load_cross_table_csv(FLCCT2PFS)

    ## load input data 
    gout = pylcu.get_grid_reg(FGOUT)
    
    ## open output dataset
    print('Opening: ',FOUTAGG)
    nc = Dataset(FOUTAGG,'r')
    chunks = nc.variables['classfr'].chunking()
    chunks = np.array(chunks)
    classfr = nc.variables['classfr']
    norm_cci = False
  else:
    ## not reading from input, need to define chunking 
    if gout['nlon']==0: 
      chunks=[1,1,np.minimum(gout['nptot'],chunkS**2)]
    else:
      chunks=[1,1,gout['nlon']]
    norm_cci = True  
    
  ## Generate the output file
  ncO = gen_nc_ifs(FOUTIFS,gout,lcI)
  
  if gout['nlon']==0:
    nlat=1
    nlon=gout['nptot']
  else:
    nlat=gout['nlat']
    nlon=gout['nlon']
    
  nslot=len(np.arange(0,nlat,chunks[1]))*len(np.arange(0,nlon,chunks[2]))
  print('LIFS chunks:',chunks,nslot) 
  ik=0
  t = time.time()
  for ilat in np.arange(0,nlat,chunks[1]):
    for ilon in np.arange(0,nlon,chunks[2]):
      if norm_cci : 
        classfrC = cci_count[:,ilat:ilat+chunks[1],ilon:ilon+chunks[2]] / cci_Ncount[ilat:ilat+chunks[1],ilon:ilon+chunks[2]]
      else:
        # reading from file 
        classfrC =  nc.variables['classfr'][:,ilat:ilat+chunks[1],ilon:ilon+chunks[2]].filled()
      ## Apply cross walking table 
      vtype=np.swapaxes(np.matmul(np.swapaxes(classfrC,0,-1),lcI['cwt']),0,-1)
      ## derived veg. cover for IFS
      ## cci derived land-sea mask, 1 - ocean
      lsm = 1.-(vtype[lcI['ifs_ind']==15,:,:]+vtype[lcI['ifs_ind']==14,:,:]).squeeze()
      land = lsm>fmin
      if np.sum(land) > 0 :
        ## Low and high vegetation cover  
        cvl = np.where(land,
                      np.sum(vtype[lcI['ifs_vtype']=='L',:,:],axis=0)/lsm,0.)
        cvh = np.where(land,
                      np.sum(vtype[lcI['ifs_vtype']=='H',:,:],axis=0)/lsm,0.)
        
        ## dominant types 
        xtmp=vtype.copy();xtmp[lcI['ifs_vtype']!='L',:,:]=-1 # set cover all other types execept low to -1                                 
        tvl = np.where(cvl>fmin,
                      lcI['ifs_ind'][np.argmax(xtmp,axis=0)],0)

        xtmp=vtype.copy();xtmp[lcI['ifs_vtype']!='H',:,:]=-1 # set cover all other types execept high to -1                                
        tvh = np.where(cvh>fmin,
                      lcI['ifs_ind'][np.argmax(xtmp,axis=0)],0)
      else:
        cvl=land*0
        cvh=cvl
        tvl=cvl
        tvh=cvl

      ## Save to the output 
      ncO.variables['lsm'][ilat:ilat+chunks[1],ilon:ilon+chunks[2]] = lsm
      ncO.variables['cvl'][ilat:ilat+chunks[1],ilon:ilon+chunks[2]] = cvl
      ncO.variables['cvh'][ilat:ilat+chunks[1],ilon:ilon+chunks[2]] = cvh
      ncO.variables['tvl'][ilat:ilat+chunks[1],ilon:ilon+chunks[2]] = tvl
      ncO.variables['tvh'][ilat:ilat+chunks[1],ilon:ilon+chunks[2]] = tvh

      ik=ik+1
      if ik%100 == 0:
          sys.stdout.write("\rProcessing :%i [%3i %%] %.1f sec"%(ik,100.*ik/nslot,time.time()-t)) ;sys.stdout.flush()
          t=time.time()   
  
  #close files
  ncO.close()
  # try to close input file in case it was opened 
  try:
    nc.close()
  except:
    pass
  
  
  
  if gout['nlon']==0:
    import eccodes as ec
    print("Writing grib files ")
    nc = Dataset(FOUTIFS,'r')
    
    ## Gaussian reduced: write grib files 
    ftemplate = open(FGOUT)
    gid_template = ec.codes_grib_new_from_file(ftemplate)
    ## ftemplate must contain LSM
    sn = ec.codes_get(gid_template,'shortName')
    if sn != 'lsm':
      raise ValueError('grib template provide is not LSM!')
    ilsm = ec.codes_get_values(gid_template)
    inotland = ilsm <= 0.5
    for cv in ['cvl','cvh','tvl','tvh']:
      foutput="%s/%s"%(FOUTP,cv)
      foutgrb= open(foutput,'wb')
      xdata = nc.variables[cv][:].squeeze()
      xdata[inotland]= 0 # set to zero if not land : lsm>0.5
      pylcu.save2grb(xdata,gid_template,foutgrb,{'shortName':cv})
      foutgrb.close()
    ec.codes_release(gid_template)
    ftemplate.close()
    nc.close()
  
  print(' ')
  print("Finished LIFS in %6.2f seconds"%(time.time()-tstart))
    
print("Finished ALL in %6.2f seconds"%(time.time()-tstart0))
