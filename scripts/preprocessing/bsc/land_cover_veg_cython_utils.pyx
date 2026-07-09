#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=False
## To compile: 
## python setup_cython.py build_ext --inplace

from __future__ import print_function
import numpy as np
cimport numpy as np



def agg_cgls_red(np.ndarray[np.float64_t,ndim=1] lat_hre,
               np.ndarray[np.float64_t,ndim=1] lon_hre,
               np.ndarray[np.float32_t,ndim=2] data_hre,
               np.ndarray[np.float32_t,ndim=2] lsm_hre,
               np.ndarray[np.float32_t,ndim=2] cvl_hre,
               np.ndarray[np.float32_t,ndim=2] cvh_hre,
               np.ndarray[np.float32_t,ndim=3] data_out,
               np.ndarray[np.float32_t,ndim=1] latr_gg,
               np.ndarray[np.int32_t,ndim=1] nplon_gg,
               np.ndarray[np.int32_t,ndim=1] nplon_ggC,
               float zmiss
               ):
  
  cdef:
    int ix_hre,iy_hre
    int ilat_gg,ilon_gg,igp_gg,npgg,ix,iy
    int nx_hre,ny_hre
    float dlon_gg,
    float latMX,latMN,lonMX,lonMN,nn,igp,nnlon
 
  nx_hre=lon_hre.shape[0]
  ny_hre=lat_hre.shape[0] 
  npgg = data_out.shape[2] 

  cdef np.ndarray[np.float32_t, ndim = 3] count = np.zeros((3,1,npgg), dtype=np.float32)
  
  nn=0
  latMX=-91
  latMN=91
  lonMX=-181
  lonMN=181
  nnlon=0
  igp=86904
  
  for iy_hre in range(ny_hre):
    # find near latitude
    ilat_gg = np.abs(latr_gg  - lat_hre[iy_hre]).argmin()
    # lon increments on that latitude 
    dlon_gg = (360./nplon_gg[ilat_gg])
    
    for ix_hre in range(nx_hre):
      # find near longitude indes
      #ilon_gg = int ( (lon_hre[ix_hre] - (0.5*dlon_gg)) / dlon_gg )
      ilon_gg = int(( 0.5*dlon_gg+lon_hre[ix_hre] )/dlon_gg)
      if (ilon_gg == nplon_gg[ilat_gg]):
        ilon_gg=0 # go to longitude zero
      if ( ilon_gg < 0 or ilon_gg >= nplon_gg[ilat_gg]):
        print(lon_hre[ix_hre],dlon_gg)
        raise ValueError('**ERROR** ilon_gg out of bounds ilon_gg=%i,max=%i'%(ilon_gg,nplon_gg[ilat_gg]-1))
      # index on gaussian reduced 
      igp_gg = nplon_ggC[ilat_gg] - nplon_gg[ilat_gg] + ilon_gg
      if ( igp_gg < 0 or igp_gg > npgg-1):
        print(ilat_gg,nplon_gg[ilat_gg],npgg,ilon_gg,igp_gg,lat_hre[iy_hre],lon_hre[ix_hre])
        raise ValueError('**ERROR** igp_gg out of bounds')
      ix=igp_gg
      iy=0
      
      if ix == igp :
        nn=nn+1
        nnlon=nplon_gg[ilat_gg]
        if lat_hre[iy_hre] > latMX: latMX=lat_hre[iy_hre]
        if lat_hre[iy_hre] < latMN: latMN=lat_hre[iy_hre]
        if lon_hre[ix_hre] > lonMX: lonMX=lon_hre[ix_hre]
        if lon_hre[ix_hre] < lonMN: lonMN=lon_hre[ix_hre]
        
     
      if data_hre[iy_hre,ix_hre] != zmiss and lsm_hre[iy_hre,ix_hre]>0.5:
        if cvl_hre[iy_hre,ix_hre] > cvh_hre[iy_hre,ix_hre]:
          data_out[0,iy,ix] = data_out[0,iy,ix] + cvl_hre[iy_hre,ix_hre]*data_hre[iy_hre,ix_hre]
          count[0,iy,ix] = count[0,iy,ix]+cvl_hre[iy_hre,ix_hre]
        else:
          data_out[1,iy,ix] = data_out[1,iy,ix] + cvh_hre[iy_hre,ix_hre]*data_hre[iy_hre,ix_hre]
          count[1,iy,ix] = count[1,iy,ix]+cvh_hre[iy_hre,ix_hre]
        data_out[2,iy,ix] = data_out[2,iy,ix] + data_hre[iy_hre,ix_hre]
        count[2,iy,ix] = count[2,iy,ix]+1
  
  print('dbg:',igp,nnlon,nn,latMN,latMX,lonMN,lonMX)    
  ## Normalize 
  for iy in range(1):
    for ix in range(npgg):
      for it in range(3):
        if count[it,iy,ix]>0  :
          data_out[it,iy,ix] = data_out[it,iy,ix] / count[it,iy,ix]
        else:
          data_out[it,iy,ix] = zmiss
      if count[0,iy,ix] == 0 and count[2,iy,ix] > 0:
        data_out[0,iy,ix] = data_out[2,iy,ix]
      if count[1,iy,ix] == 0 and count[2,iy,ix] > 0:
        data_out[1,iy,ix] = data_out[2,iy,ix]
       
def agg_cgls_reg(np.ndarray[np.float64_t,ndim=1] lat_hre,
               np.ndarray[np.float64_t,ndim=1] lon_hre,
               np.ndarray[np.float32_t,ndim=2] data_hre,
               np.ndarray[np.float32_t,ndim=2] lsm_hre,
               np.ndarray[np.float32_t,ndim=2] cvl_hre,
               np.ndarray[np.float32_t,ndim=2] cvh_hre,
               np.ndarray[np.float32_t,ndim=3] data_out,
               double lat0_out,
               double dlat_out,
               double lon0_out,
               double dlon_out,
               int nlat_out,
               int nlon_out,
               int subarea,
               float zmiss
               ):
  
  cdef:
    int ix_hre,iy_hre,ix,iy
    int nx_hre,ny_hre
   
  nx_hre=lon_hre.shape[0]
  ny_hre=lat_hre.shape[0]  

  cdef np.ndarray[np.int32_t, ndim = 1] ixA = np.zeros((nx_hre), dtype=np.int32)
  cdef np.ndarray[np.float32_t, ndim = 3] count = np.zeros((3,nlat_out,nlon_out), dtype=np.float32)
  
  
  ## precompute ix as it is not necessary to compute it for each latitude
  for ix_hre in range(nx_hre):
      # near longitude 
      ix = int( ( lon_hre[ix_hre] - (lon0_out-0.5*dlon_out) ) / dlon_out ) 
      ixA[ix_hre] = ix
      if ( ix < 0 or ix > (nlon_out -1 ) ):
        if subarea == 1:
          ixA[ix_hre] = -1
          continue
        else:
          raise ValueError('**ERROR** ix out of bounds indlon_hre[ix_hre]=[%i],%i'%(ix_hre,ix))
  
  # do the count 
  for iy_hre in range(ny_hre):
    # near latitude 
    iy = int( ( lat_hre[iy_hre] - (lat0_out-0.5*dlat_out) ) / dlat_out )
    if ( iy < 0 or iy > (nlat_out -1 ) ):
      if subarea == 1 :
        continue
      else:
        raise ValueError('**ERROR** iy out of bounds=[%i],%i'%(iy_hre,iy))
    
    for ix_hre in range(nx_hre):
      # near longitude 
      ix = ixA[ix_hre]
      if ix == -1 : 
        continue
      
      if data_hre[iy_hre,ix_hre] != zmiss and lsm_hre[iy_hre,ix_hre]>0.5:
        if cvl_hre[iy_hre,ix_hre] > cvh_hre[iy_hre,ix_hre]:
          data_out[0,iy,ix] = data_out[0,iy,ix] + cvl_hre[iy_hre,ix_hre]*data_hre[iy_hre,ix_hre]
          count[0,iy,ix] = count[0,iy,ix]+cvl_hre[iy_hre,ix_hre]
        else:
          data_out[1,iy,ix] = data_out[1,iy,ix] + cvh_hre[iy_hre,ix_hre]*data_hre[iy_hre,ix_hre]
          count[1,iy,ix] = count[1,iy,ix]+cvh_hre[iy_hre,ix_hre]
        data_out[2,iy,ix] = data_out[2,iy,ix] + data_hre[iy_hre,ix_hre]
        count[2,iy,ix] = count[2,iy,ix]+1
      
  ## Normalize 
  for iy in range(nlat_out):
    for ix in range(nlon_out):
      for it in range(3):
        if count[it,iy,ix]>0  :
          data_out[it,iy,ix] = data_out[it,iy,ix] / count[it,iy,ix]
        else:
          data_out[it,iy,ix] = zmiss
      if count[0,iy,ix] == 0 and count[2,iy,ix] > 0:
        data_out[0,iy,ix] = data_out[2,iy,ix]
      if count[1,iy,ix] == 0 and count[2,iy,ix] > 0:
        data_out[1,iy,ix] = data_out[2,iy,ix]
  
  
def agg_cci_reg(
               np.ndarray[np.float64_t,ndim=1] lat_hre,
               np.ndarray[np.float64_t,ndim=1] lon_hre,
               np.ndarray[np.uint8_t,ndim=2] data_hre,
               np.ndarray[np.uint16_t,ndim=3] data_out,
               np.ndarray[np.uint16_t,ndim=1] lcci_ind,
               float lat0_out,
               float dlat_out,
               float lon0_out,
               float dlon_out,
               int nlat_out,
               int nlon_out,
               int subarea,
               ):
  cdef:
    int ix_hre,iy_hre,ix,iy
    int nx_hre,ny_hre,ind,nclass_max
  nx_hre=lon_hre.shape[0]
  ny_hre=lat_hre.shape[0]
  nclass_max=data_out.shape[0]
  
  cdef np.ndarray[np.int32_t, ndim = 1] ixA = np.zeros((nx_hre), dtype=np.int32) 
  
  ## precompute ix as it is not necessary to compute it for each latitude
  for ix_hre in range(nx_hre):
      # near longitude 
      ix = int( ( lon_hre[ix_hre] - (lon0_out-0.5*dlon_out) ) / dlon_out ) 
      ixA[ix_hre] = ix
      if ( ix < 0 or ix > (nlon_out -1 ) ):
        if subarea == 1:
          ixA[ix_hre] = -1
          continue
        else:
          raise ValueError('**ERROR** ix out of bounds indlon_hre[ix_hre]=[%i],%i'%(ix_hre,ix))
      
  # do the count 
  for iy_hre in range(ny_hre):
    # near latitude 
    iy = int( ( lat_hre[iy_hre] - (lat0_out-0.5*dlat_out) ) / dlat_out )
    if ( iy < 0 or iy > (nlat_out -1 ) ):
      if subarea == 1 :
        continue
      else:
        raise ValueError('**ERROR** iy out of bounds=[%i],%i'%(iy_hre,iy))
    
    for ix_hre in range(nx_hre):
      # near longitude 
      ix = ixA[ix_hre]
      if ix == -1 : 
        continue
      # esa cci land cover index 
      ind=lcci_ind[data_hre[iy_hre,ix_hre]]
      if (ind < 0 or ind > nclass_max ):
        raise ValueError('Ind out of bounds',iy_hre,ix_hre,iy,ix,data_hre[iy,ix],ind)
      # accumulate
      #if (np.abs(lat_hre[iy_hre]-40)<0.01 and np.abs(lon_hre[ix_hre]--5)<0.01):
        #print(iy_hre,ix_hre,iy,ix,data_hre[iy,ix],ind)
      data_out[ind,iy,ix] = data_out[ind,iy,ix] + 1 
      
def agg_cci_red(
               np.ndarray[np.float64_t,ndim=1] lat_hre,
               np.ndarray[np.float64_t,ndim=1] lon_hre,
               np.ndarray[np.uint8_t,ndim=2] data_hre,
               np.ndarray[np.uint16_t,ndim=3] data_out,
               np.ndarray[np.uint16_t,ndim=1] lcci_ind,
               np.ndarray[np.float32_t,ndim=1] latr_gg,
               np.ndarray[np.int32_t,ndim=1] nplon_gg,
               np.ndarray[np.int32_t,ndim=1] nplon_ggC):
  cdef:
    int ix_hre,iy_hre
    int ilat_gg,ilon_gg,igp_gg,npgg
    int nx_hre,ny_hre,nclass_max,ind,
    float dlon_gg
    
  nx_hre=lon_hre.shape[0]
  ny_hre=lat_hre.shape[0]
  npgg = data_out.shape[2]
  nclass_max=data_out.shape[0]
  
  for iy_hre in range(ny_hre):
    # find near latitude
    ilat_gg = np.abs(latr_gg  - lat_hre[iy_hre]).argmin()
    # lon increments on that latitude 
    dlon_gg = (360./nplon_gg[ilat_gg])
    
    for ix_hre in range(nx_hre):
      # find near longitude indes
      #ilon_gg = int ( (lon_hre[ix_hre] - (0.5*dlon_gg)) / dlon_gg )
      ilon_gg = int(( 0.5*dlon_gg+lon_hre[ix_hre] )/dlon_gg)
      if (ilon_gg == nplon_gg[ilat_gg]):
        ilon_gg=0 # go to longitude zero
      if ( ilon_gg < 0 or ilon_gg >= nplon_gg[ilat_gg]):
        print(lon_hre[ix_hre],dlon_gg)
        raise ValueError('**ERROR** ilon_gg out of bounds ilon_gg=%i,max=%i'%(ilon_gg,nplon_gg[ilat_gg]-1))
      # index on gaussian reduced 
      igp_gg = nplon_ggC[ilat_gg] - nplon_gg[ilat_gg] + ilon_gg
      if ( igp_gg < 0 or igp_gg > npgg-1):
        print(ilat_gg,nplon_gg[ilat_gg],npgg,ilon_gg,igp_gg,lat_hre[iy_hre],lon_hre[ix_hre])
        raise ValueError('**ERROR** igp_gg out of bounds')
      
      ind=lcci_ind[data_hre[iy_hre,ix_hre]]
      if (ind < 0 or ind > nclass_max ):
        raise ValueError('Ind out of bounds',iy_hre,ix_hre,data_hre[iy_hre,ix_hre],ind)
      ## accumulate 
      data_out[ind,0,igp_gg] = data_out[ind,0,igp_gg] + 1 
      
 
