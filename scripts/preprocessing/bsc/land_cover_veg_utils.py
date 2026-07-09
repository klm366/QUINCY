from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import pandas as pd 

try:
  import eccodes as ec
except:
  print("eccodes not loaded...")

PLANET_RADIUS=6371000.  # default radius for area calculations


## function to save fields to grib output 
def save2grb(values,gid,fout,pparams={}):
  print('writing to:',fout)
  gidC = ec.codes_clone(gid)
  for kk in pparams.keys():
    ec.codes_set(gidC,kk,pparams[kk])
  ec.codes_set_values(gidC,values)
  if type(fout) is str:
    f = open(fname,'wb')
    ec.codes_write(gidC,fout)
    ec.codes_release(gidC)
    f.close()
  else:
    ec.codes_write(gidC,fout)
    ec.codes_release(gidC)

def get_grid_reg(fin,verbose=True,retbounds=False,retgarea=False):
  PLANET_RADIUS=6371000.  # default radius for area calculations
  """ 
  Get grib grid information 
  grid_info = get_grib_grid(fin,verbose)
  fin - netcdf file (will read the first field only)
  verbose - if true prints some information
  grid_info: dictionary with the output:
    grid_info['nptot'] = npG # total number of grid points 
    grid_info['nlat'] = nlatG # number of latitude lines
    grid_info['nlon'] = number of longitude lines 
    grid_info['latr'] = latrG  # gaussian regular latitude
    grid_info['lonr'] = lonrG  # gaussian regular longitude
    grid_info['dlon'] = dlon  # longitude increment in each latitude line
    grid_info['weights'] =   # gaussian weigths
  """
  if verbose:
    print('get_grid_reg: reading file',fin)
  isNC=True
  try:
    nc =  Dataset(fin,'r')
  except:
    isNC=False
  if isNC: 
    ## information from grib file
    nlat=None
    for cdim in ['lat','latitude']:
      if cdim in nc.dimensions:
        nlat = len(nc.dimensions[cdim])
        latr = nc.variables[cdim][:]
        break
    nlon=None
    for cdim in ['lon','longitude']:
      if cdim in nc.dimensions:
        nlon = len(nc.dimensions[cdim])
        lonr = nc.variables[cdim][:]
        break
    if nlat is None or nlon is None:
      raise ValueError('get_grid_reg: Could not find latitude or longitude nlat,nlon',nlat,nlon)
    nc.close()
  else:
    f = open(fin,'rb')
    gid = ec.codes_grib_new_from_file(f)
    gtype = ec.codes_get(gid,'gridType')
    if gtype in ["reduced_gg","regular_gg"]:
      ec.codes_release(gid)
      f.close()
      grib_info = get_grib_grid(fin,verbose=verbose,retbounds=retbounds)
      return grib_info
      
    nlat = ec.codes_get(gid,'Nj')
    nlon = ec.codes_get(gid,'Ni')
    lat0 = ec.codes_get(gid,'latitudeOfFirstGridPointInDegrees')
    lon0 = ec.codes_get(gid,'longitudeOfFirstGridPointInDegrees')
    lat1 = ec.codes_get(gid,'latitudeOfLastGridPointInDegrees')
    lon1 = ec.codes_get(gid,'longitudeOfLastGridPointInDegrees')
    latr = np.linspace(lat0,lat1,nlat)
    lonr = np.linspace(lon0,lon1,nlon)
    
    f.close()
    
  if len(lonr)>1:
    dlonG = ( lonr[-1]-lonr[0])/ float(nlon-1) 
  else:
    dlonG=1
  if len(latr)>1:
    dlatG = ( latr[-1]-latr[0])/ float(nlat-1) 
  else:
    dlatG=1

  grid_info = {}
  grid_info['nptot'] = nlat*nlon # total number of grid points 
  grid_info['nlat'] = nlat # number of latitude lines
  grid_info['nlon'] = nlon # to identify we'r gaussian !
  
  grid_info['latr'] = latr.astype('f8')  # regular latitude
  grid_info['lonr'] = lonr.astype('f8')  # regular longitude
  grid_info['dlon'] = dlonG  # longitude increments 
  grid_info['dlat'] = dlatG  # latitude increments 
  
  
  if retbounds:
    xbounds=np.zeros((npG,4),dtype='f4')  # west,east,east,west
    ybounds=np.zeros((npG,4),dtype='f4')  # south,south,north,north
    #
    #       3 -----  2
    #         \   \
    #         \   \
    #       0 -----  1
    #

    ybounds[:,0]=latG+dlatG*.5
    ybounds[:,1]=latG+dlatG*.5
    ybounds[:,2]=latG-dlatG*.5
    ybounds[:,3]=latG-dlatG*.5
    xbounds[:,0]=lonG-dlonG*.5
    xbounds[:,1]=lonG+dlonG*.5
    xbounds[:,2]=lonG+dlonG*.5
    xbounds[:,3]=lonG-dlonG*.5
    

    grid_info['lon_bnds'] = xbounds  # longitude bounds
    grid_info['lat_bnds'] = ybounds  # longitude bounds

    ## area computations
  if retgarea:
    PLANET_RADIUS=6371000.
    area= np.empty((nlat),dtype='f8')
    ww=np.cos(latr*np.pi/180.)
    ww=ww/np.sum(ww)
    for ilat in xrange(nlat):
      area[ilat] = 4*np.pi*PLANET_RADIUS**2*ww[ilat]/nlon
    grid_info['area'] = area
  

  if verbose:
    for kk in ['nlat','dlat','nlon','dlon']:
      print(kk,grid_info[kk])
    print('lat0,lat1',grid_info['latr'][0],grid_info['latr'][-1])
    print('lon0,lon1',grid_info['lonr'][0],grid_info['lonr'][-1])
  return grid_info

def load_cross_table_csv(finput):
  ## load data as pandas data frame, assuming it's proper csv ... 
  df = pd.read_csv(finput)
  
  ## empty dictionary to store information 
  cwtI={}

  ## CCI indexes in uint8 
  cwtI['cci_ind'] = df['CCI_IND'][2:].to_numpy('uint8')
  ## CCI land cover names
  cwtI['cci_name'] = df['CCI_NAME'][2:].to_list()
  ## IFS vegetation types name 
  cwtI['ifs_name'] = df.columns[3:].to_list()
  ## IFS vegetation index 
  cwtI['ifs_ind'] = df.iloc[0][3:].to_numpy('i4')  
  ## IFS veg. type H, L, N.A 
  cwtI['ifs_vtype'] = df.iloc[1][3:].to_numpy()
  ## number of CCI classes
  cwtI['ncci'] = len(cwtI['cci_ind'])
  ## number veg. types IFS
  cwtI['nifs'] = len(cwtI['ifs_ind'])
  ## Cross walking table 
  cwtI['cwt'] = df.iloc[2:,3:].fillna(0).to_numpy('f4')/100. 
  ## check that cwt sum 100%
  for i in range(cwtI['ncci']):
    xsum = np.sum(cwtI['cwt'][i,:])
    if xsum != 1: 
      print("%i %s weights do not sum to 1: %f"%(
        cwtI['cci_ind'][i],cwtI['cci_name'][i],xsum))
  return cwtI

def find_bbox(bbox,ginfo,lat=None,lon=None,verbose=False):
  bboxO=[-999,-999,-999,-999]

  il0=np.int((bbox[0]-ginfo['lon0'])/ginfo['dlon'])
  il1=np.int((bbox[1]-ginfo['lon0'])/ginfo['dlon'])
  if (il1>il0):
    bboxO[0]=il0
    bboxO[1]=il1
  else:
    bboxO[1]=il0
    bboxO[0]=il1
  
  if verbose:
    il0=bboxO[0]
    il1=bboxO[1]
    print('lon ',il0,il1,lon[il0],lon[il1],bbox[0],bbox[1],len(lon[il0:il1+1]))
  il0=np.int((bbox[2]-ginfo['lat0'])/ginfo['dlat'])
  il1=np.int((bbox[3]-ginfo['lat0'])/ginfo['dlat'])
  if (il1>il0):
    bboxO[2]=il0
    bboxO[3]=il1
  else:
    bboxO[3]=il0
    bboxO[2]=il1
  if verbose:  
    il0=bboxO[2]
    il1=bboxO[3]
    print('lat ',il0,il1,lat[il0],lat[il1],bbox[2],bbox[3],len(lat[il0:il1+1]))
  return bboxO

def get_grib_grid(fin,verbose=True,retbounds=False):
  """ 
  Get grib grid information 
  grid_info = get_grib_grid(fin,verbose)
  fin - grib file (will read the first field only)
  verbose - if true prints some information
  grid_info: dictionary with the output:
    grid_info['nptot'] = npG # total number of grid points 
    grid_info['nlat'] = nlatG # number of latitude lines
    grid_info['nlon'] = 0 # to identify we'r gaussian !
    grid_info['nplon'] = nplonG # number longitudes in each latitude 
    grid_info['lat'] = latG # global lat
    grid_info['lon'] = lonG # global lon
    grid_info['nplonC'] = nplonGC  # cumulative number of points in each latitude line 
    grid_info['latr'] = latrG  # gaussian regular latitude
    grid_info['lonr'] = lonrG  # gaussian regular longitude
    grid_info['dlon'] = dlon  # longitude increment in each latitude line
    grid_info['weights'] =   # gaussian weigths
  """
  if verbose:
    print('get_grib_grid: reading file',fin)
  f = open(fin,'rb')
  gid = ec.codes_grib_new_from_file(f)
  
  ## information from grib file 
  nlatG = ec.codes_get(gid,'Nj') # number of latitudes
  gtype = ec.codes_get(gid,'gridType')
  if gtype == "reduced_gg":
    nplonG = ec.codes_get_array(gid,'pl')  # number of point in each latitude
  elif gtype == "regular_gg":
    nlon = ec.codes_get_array(gid,'Ni')
    nplonG = np.zeros((nlatG),dtype='i4')+nlon  
  
  ec.codes_release(gid)
  f.close()
  ## derived
  nplonGC = np.cumsum(nplonG,dtype='i4')
  npG = np.sum(nplonG,dtype='i4')
  lonr = np.linspace(0.,360.,np.max(nplonG),endpoint=False)
  dlonG = 360./nplonG
  ## gen gaussian lat and weights
  latr,weights = get_reduced_lats_weights(nlatG)
  
  latG = np.empty(npG,dtype='f8')
  lonG = np.empty(npG,dtype='f8')
  istart=0
  for ilat in range(nlatG):
    iend = istart+nplonG[ilat]
    latG[istart:iend] = latr[ilat]
    lonG[istart:iend] = np.linspace(0.,360.,nplonG[ilat],endpoint=False)
    istart = iend


  grid_info = {}
  grid_info['nptot'] = npG # total number of grid points 
  grid_info['nlat'] = nlatG # number of latitude lines
  grid_info['nlon'] = 0 # to identify we'r gaussian !
  grid_info['nplon'] = nplonG # number longitudes in each latitude 
  grid_info['lat'] = latG # global lat
  grid_info['lon'] = lonG # global lon
  grid_info['nplonC'] = nplonGC  # cumulative number of points in each latitude line 
  grid_info['latr'] = latr  # gaussian regular latitude
  grid_info['lonr'] = lonr  # gaussian regular longitude
  grid_info['dlon'] = dlonG  # longitude increment in each latitude line
  grid_info['gweights'] = weights  # longitude increment in each latitude line
  
  
  if retbounds:
    xbounds=np.zeros((npG,4),dtype='f4')  # west,east,east,west
    ybounds=np.zeros((npG,4),dtype='f4')  # south,south,north,north
    #
    #       3 -----  2
    #         \   \
    #         \   \
    #       0 -----  1
    #

    ## compute lat edges
    latedges = np.zeros((nlatG+1),dtype='f8')
    latedges[0] = 90.
    latedges[1:-1] = np.arcsin(1-np.cumsum(weights[:-1]))*180/np.pi
    latedges[-1] = -90.
    ## filll the arrays 
    istart=0
    for ilat in xrange(nlatG):
      iend = istart+nplonG[ilat]
      ybounds[istart:iend,0]=latedges[ilat+1]
      ybounds[istart:iend,1]=latedges[ilat+1]
      ybounds[istart:iend,2]=latedges[ilat]
      ybounds[istart:iend,3]=latedges[ilat]
      xbounds[istart:iend,0]=lonG[istart:iend]-dlonG[ilat]*.5
      xbounds[istart:iend,1]=lonG[istart:iend]+dlonG[ilat]*.5
      xbounds[istart:iend,2]=xbounds[istart:iend,1]
      xbounds[istart:iend,3]=xbounds[istart:iend,0]
      istart = iend

    grid_info['lon_bnds'] = xbounds  # longitude bounds
    grid_info['lat_bnds'] = ybounds  # longitude bounds

    ## area computations
  area= np.empty(npG,dtype='f8')
  #area = np.empty(npG,dtype='f8')
  istart=0
  for ilat in range(nlatG):
    iend = istart+nplonG[ilat]
    #xarea = 0.5*np.abs(np.sin(np.deg2rad(ybounds[istart,0]))-
                       #np.sin(np.deg2rad(ybounds[istart,2])) )/nplonG[ilat]
    xarea = 4*np.pi*PLANET_RADIUS**2*weights[ilat] / (nplonG[ilat]*np.sum(weights))
    area[istart:iend] = xarea
    #area1[istart:iend] = xarea1 
    istart = iend

  grid_info['area'] = area  # longitude increment in each latitude line
  #grid_info['area'] = area  # longitude increment in each latitude line

  if verbose:
    print("nptot,nlat,nlon:",' ',grid_info['nptot'],grid_info['nlat'],grid_info['nlon'])
    print('get_grib_grid: finished')
  return grid_info

def get_reduced_lats_weights(N,FLOC=''):
  """ 
  Compute the latitudes and weigths of a gaussian grid with N latitudes:
  lats,weights = get_reduced_lats_weights(N) 
  
  Using:
  Legendre quadrature points on [-1,1]. Returns the roots of Nth order Legendre polynomial P_N(x) and   
  weights (w) to use in Gaussian Quadrature over [-1,1] with weighting function 1.
  
  Adapted from  Anna Agusti-Panareda
  """
  #import scipy.special
  #sinlat,weight = scipy.special.orthogonal.p_roots(N)
  import os.path
  fname=FLOC+'gauss_%i.npy'%(N/2)
  if os.path.isfile(fname):
    OUT = np.load(fname)
    lat = OUT[:,0]
    weight = OUT[:,1]
  else:
    OUT=np.zeros((N,2))
    OUT[:,0],OUT[:,1] = np.polynomial.legendre.leggauss(N)
  
    lat=(np.arcsin(OUT[:,0])*180.0/np.pi)[::-1] # invert to be between 90: -90
    weight = OUT[:,1][::-1]
  return lat,weight

def wrapTo360(xin):
  """ 
  wrap longitues between 0 to 360
  xout = wrapTo360(xin)
  xin: input np.array with longitudes 
  xout: output np.array with longitudes 
  """
  xout = xin.copy()
  if np.size(xout) > 1 :
      xout[xin<0.] = 360. + xin[xin<0.]
  else:
      if xin<0. :
          xout = 360. + xin 
  return xout 
