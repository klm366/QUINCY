##    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
##    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
##    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
##
##    This file is part of the QUINCY model.
##    The QUINCY model is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License version 3,
##    as published by the Free Software Foundation.
##
##    The QUINCY model is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.

##
## script to generate surface run files for SLM
##

## library
library(ncdf4)

#########################################################################################
## run this script by:
## cd contrib/forcing_scripts/
## Rscript ./produce_face_site_info_file.R


#########################################################################################
## output 
## the output files will be written in this folder


#########################################################################################
## site_list.dat containing basic site-level data
## input paths for two files
ipath_nc         <- "/Net/Groups/BGI/model/OCN/input/point/FACE/"
ipath_info       <- "/Net/Groups/BGI/people/szaehle/Projects/QUINCY/GenerateFACEForcing_2017/"
allsiteinfo      <- read.table(paste0(ipath_info,"site_list.dat"),stringsAsFactor=F,header=T)

## get site names and start/end years of data
allsitestartyear <- allsiteinfo$Start
allsiteendyear   <- allsiteinfo$End
allsitename      <- allsiteinfo$Site
allsite_newname  <- allsiteinfo$NEWNAME   # DK-AMB DK-ELE OR-AMB OR-ELE
site_list        <- allsitename           # DUKE   DUKE   ORNL   ORNL
pft              <- allsiteinfo$PFT

## create empty vectors for output
lon_orig <- numeric(length(site_list))
lat_orig <- numeric(length(site_list))

## loop over the four sites and get their lon/lat (from the two files DUKE* & ORNL* _FACE.nc)
for (s in 1:length(site_list)) {
  fid         <- nc_open(paste0(ipath_nc,site_list[s],"_FACE.nc"))
  lon_orig[s] <- ncvar_get(fid,"x")
  lat_orig[s] <- ncvar_get(fid,"y")
  nc_close(fid)
}

## get site information and save into vector
gmt_ref   <- allsiteinfo[,15]
sand      <- allsiteinfo[,12]
silt      <- allsiteinfo[,11]
clay      <- 1 - sand - silt
awc       <- allsiteinfo[,13]
ph        <- allsiteinfo[,14]
lai       <- allsiteinfo[,10]*0 + -9999  # i.e., set to -9999
nleaf     <- allsiteinfo[,11]*0 + -9999
sla       <- allsiteinfo[,12]*0 + -9999
height    <- allsiteinfo[,13]*0 + -9999
plantyear <- allsiteinfo[,16]


#########################################################################################
## get supplementary data from Soil Grids Database
## new input path
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/"

# note:
# sand silt awc are read & used from the above site_list.dat file

ncid  <- nc_open(paste0(ipath,"SNDPPT.soilgrid.3600.1800.7.nc"))
lon_sg<-ncvar_get(ncid,'longitude')
lat_sg<-ncvar_get(ncid,'latitude')
#sand_sg<-ncvar_get(ncid,'SNDPPT')
nc_close(ncid)

#ncid<-nc_open(paste(ipath,'SLTPPT.soilgrid.3600.1800.7.nc',sep=''))
#silt_sg<-ncvar_get(ncid,'SLTPPT')
#nc_close(ncid)

ncid<-nc_open(paste(ipath,'BLDFIE.soilgrid.3600.1800.7.nc',sep=''))
bulk_density_sg<-ncvar_get(ncid,'BLDFIE')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'PHIHOX.soilgrid.3600.1800.7.nc',sep=''))
ph_sg<-ncvar_get(ncid,'PHIHOX')
nc_close(ncid)

#ncid<-nc_open(paste(ipath,'AWCtS.soilgrid.3600.1800.7.nc',sep=''))
#fc_vol_sg<-ncvar_get(ncid,'AWCtS')
#nc_close(ncid)

#ncid<-nc_open(paste(ipath,'WWP.soilgrid.3600.1800.7.nc',sep=''))
#pwp_vol_sg<-ncvar_get(ncid,'WWP')
#nc_close(ncid)

ncid<-nc_open(paste(ipath,'BDRICM.soilgrid.3600.1800.nc',sep=''))
depth_to_bedrock_sg<-ncvar_get(ncid,'BDRICM')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'TAXOUSDA_10km.soilgrid.3600.1800.nc',sep=''))
taxousda_sg<-ncvar_get(ncid,'TAXOUSDA_10km')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'TAXNWRB.soilgrid.3600.1800.nc',sep=''))
taxnwrb_sg<-ncvar_get(ncid,'TAXNWRB')
nc_close(ncid)

## create new (empty) vectors/arrays for output
depth_begin  <- c(0.0,0.05,0.15,0.3,0.6,1.0,2.0)
depth_end    <- c(0.05,0.15,0.3,0.6,1.0,2.0,100)
depth_sl     <- depth_end - depth_begin
taxnwrb      <- numeric(length(site_list))
taxousda     <- numeric(length(site_list))
depth        <- numeric(length(site_list))
bulk_density <- numeric(length(site_list))
awc_sl       <- array(NA,c(length(site_list),7))

## loop over sites: get soil data
for (s in 1:length(site_list)) {
  lonidx          <- which.min(abs(lon_sg-lon_orig[s]))
  latidx          <- which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  # depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]

  if(ph[s]==6.5) ph[s]<-mean(ph_sg[lonidx,latidx,])/10

  #  note:
  #  sand silt awc are read & used from the above site_list.dat file
  # if(is.na(sand[s])) sand[s]<-weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  # if(is.na(silt[s])) silt[s]<-weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  # awc_sl[s,]<-(fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  # if(is.na(awc[s]))  awc[s]<-weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}


#########################################################################################
## read lithology map
## new input file
ncid <- nc_open("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc")

lon_glim        <- ncvar_get(ncid,'lon')
lat_glim        <- ncvar_get(ncid,'lat')
glim_class      <- ncvar_get(ncid,'GLim')
glim_class_site <- numeric(length(site_list))
nc_close(ncid)

## loop over sites: get soil data
for (s in 1:length(site_list)) {
  lonidx             <- which.min(abs(lon_glim-lon_orig[s]))
  latidx             <- which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s] <- glim_class[lonidx,latidx]
}

## replace NAN with default values
sand[is.na(sand)] <- 0.4
silt[is.na(silt)] <- 0.4
clay              <- 1 - sand - silt
ph[is.na(ph)]     <- 6
awc[is.na(awc)]   <- 200
bulk_density[is.na(bulk_density)]       <- 1500
taxousda[is.na(taxousda)]               <- 30
taxnwrb[is.na(taxnwrb)]                 <- 27
glim_class_site[is.na(glim_class_site)] <- 2


#########################################################################################
## read site-specific data on Phosphorous
## multiple new input files
inPhosPath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/Phosphorous/v2014_06/Data/"

# soilP_depth_sg is 0.5 [m], same value for all sites

ncid             <- nc_open(paste0(inPhosPath,"Labile_Inorganic_P.nc"))
# same lon/lat as above, reading it anyway to improve modularity & reduce dependencies 
lon_sg           <- ncvar_get(ncid,'longitude')
lat_sg           <- ncvar_get(ncid,'latitude')
soilP_labile_sg  <- ncvar_get(ncid,"Labile_Inorganic_P")
nc_close(ncid)

ncid             <- nc_open(paste0(inPhosPath,"Seconday_Mineral__P.nc"))
soilP_slow_sg    <- ncvar_get(ncid,"Seconday_Mineral__P")
nc_close(ncid)

ncid             <- nc_open(paste0(inPhosPath,"Occluded_P.nc"))
soilP_occlud_sg  <- ncvar_get(ncid,"Occluded_P")
nc_close(ncid)

ncid             <- nc_open(paste0(inPhosPath,"Apatite_P.nc"))
soilP_primary_sg <- ncvar_get(ncid,"Apatite_P")
nc_close(ncid)

## create empty vectors for site-data
soilP_depth_site   <- numeric(length(site_list))
soilP_labile_site  <- numeric(length(site_list))
soilP_slow_site    <- numeric(length(site_list))
soilP_occlud_site  <- numeric(length(site_list))
soilP_primary_site <- numeric(length(site_list))

## loop over sites: get soil P data from nearest location
for (s in 1:length(site_list)) {
  # get neares lon/lat in soil map
  lonidx                <- which.min(abs(lon_sg-lon_orig[s]))
  latidx                <- which.min(abs(lat_sg-lat_orig[s]))

  soilP_depth_site[s]   <- 0.5
  soilP_labile_site[s]  <- round(soilP_labile_sg[lonidx,latidx],4)
  soilP_slow_site[s]    <- round(soilP_slow_sg[lonidx,latidx],4)
  soilP_occlud_site[s]  <- round(soilP_occlud_sg[lonidx,latidx],4)
  soilP_primary_site[s] <- round(soilP_primary_sg[lonidx,latidx],4)

  # debug output
  #print(paste0("site: ",s," ",soilP_labile_site[s]," ",soilP_slow_site[s]," ",soilP_occlud_site[s]," ",soilP_primary_site[s]))
}


#########################################################################################
## get site-specific data on qmax_org_fine_particle
## qmax_org_fine_particle was calculated by Lin (lyu@bgc-jena.mpg.de) per NWRB category available from TAXNWRB.soilgrid.3600.1800.nc
## two new input files
inNWRBcatFile     <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/TAXNWRB.soilgrid.3600.1800.nc"
inCalcQmaxorgFile <- "qmax_org_values_per_nwrb_category_20180515.csv"


ncid         <- nc_open(inNWRBcatFile)  # grid with 3600x1800 resolution
lon_sg       <- ncvar_get(ncid,'longitude')
lat_sg       <- ncvar_get(ncid,'latitude')
nwrbCat_sg   <- ncvar_get(ncid,"TAXNWRB")
nc_close(ncid)

qmaxorgTable <- read.table(inCalcQmaxorgFile, header=T, sep="", stringsAsFactors=F)

## create empty vectors for site-data
qmaxorg_site <- numeric(length(site_list))

## loop over sites: get NWRB cat data from nearest location & get according Qmax_org
for (s in 1:length(site_list)) {
  # get neares lon/lat in soil map
  lonidx <- which.min(abs(lon_sg-lon_orig[s]))
  latidx <- which.min(abs(lat_sg-lat_orig[s]))
  qmaxorg_site[s] <- qmaxorgTable$qmax_org_value[ as.integer(nwrbCat_sg[lonidx,latidx]) ]

  # round to 4 decimal places
  qmaxorg_site[s] <- round(qmaxorg_site[s],4)
  
  # debug output
  #print(paste0("site: ",s," ",qmaxorg_site[s]))
}



#########################################################################################
## output
## create empty output array and paste data into it
siteParam <- 27
siteinfo  <- array(NA,c(length(site_list),siteParam))

siteinfo[,1]<-allsite_newname
siteinfo[,2]<-lon_orig
siteinfo[,3]<-lat_orig
siteinfo[,4]<-pft
siteinfo[,5]<-allsitestartyear
siteinfo[,6]<-allsiteendyear
siteinfo[,7]<-gmt_ref
siteinfo[,8]<-clay
siteinfo[,9]<-silt
siteinfo[,10]<-sand
siteinfo[,11]<-awc
siteinfo[,12]<-bulk_density
siteinfo[,13]<-ph
siteinfo[,14]<-taxousda
siteinfo[,15]<-taxnwrb
siteinfo[,16]<-glim_class_site
siteinfo[,17]<-lai
siteinfo[,18]<-nleaf
siteinfo[,19]<-sla
siteinfo[,20]<-height
siteinfo[,21]<-plantyear
siteinfo[,22]<-soilP_depth_site
siteinfo[,23]<-soilP_labile_site
siteinfo[,24]<-soilP_slow_site
siteinfo[,25]<-soilP_occlud_site
siteinfo[,26]<-soilP_primary_site
siteinfo[,27]<-qmaxorg_site


## output
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')

# 1996-2008  (DUKE: 1996-2007 & ORNL: 1998-2008)
write.table(siteinfo,file=paste0("face_sites_1996-2008_list",Sys.Date(),".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1901
siteinfo[,6]<-2015
# 1901-2015
write.table(siteinfo,file=paste0("face_sites_1901-2015_list",Sys.Date(),".dat"),quote=F,row.names=F,col.names=T)


print(paste0("Script has finished. The new face_sites_list.dat with ",siteParam," Parameters were created."))



## end ##
