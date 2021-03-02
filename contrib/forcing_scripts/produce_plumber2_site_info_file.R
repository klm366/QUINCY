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
## script to generate site-info (a.k.a. surface or ic/bc) files for QUINCY - PLUMBER2 site-set
## PLUMBER is a model intercomparison project - the organisers have collected data from multiple sources, as indicated by the 'origin' information
## here we use the 2nd version and site-set of this project (https://modelevaluation.org)
##

## library
library(ncdf4)

#########################################################################################
## run this script by:
##        cd contrib/forcing_scripts/
##        Rscript --vanilla ./produce_plumber2_site_info_file.R


#########################################################################################
## output 
## the output files will be written into this folder


#########################################################################################
## *site_list.dat containing basic site-level data
# input path for met-forcing netcdf files
ipath            <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/PLUMBER2_siteset_data/"
## input file pattern definitions
infile_ext       <- '_Met.nc'
# get a list of all files with original met data from the ipath directory
sitesetFiles     <- Sys.glob(paste0(ipath, "*", infile_ext))
# site list info file (modified version of the file created by Silvia Caldararu <scaldra@bgc-jena.mpg.de> December 2020)
allsiteinfo      <- read.table("generate_plumber2_site_list.dat", header=T, sep="", stringsAsFactor=F)

## get site names and start/end years of data
allsitestartyear <- allsiteinfo$start
allsiteendyear   <- allsiteinfo$end
allsitename      <- allsiteinfo$site_id
allsiteorigin    <- allsiteinfo$origin      # refers to the flux tower network and version that the data is from
pft              <- allsiteinfo$PFT
site_list        <- allsitename             # "AR-SLu" "AT-Neu" ... "ZA-Kru" "ZM-Mon"

## set default NA-value for variables no data is available (for all sites)
lai       <- -9999
nleaf     <- -9999
sla       <- -9999
height    <- -9999
plantyear <- 1500 

## create empty vectors for output
lon_orig <- numeric(length(site_list))
lat_orig <- numeric(length(site_list))

# get lon/lat from the metforcing.nc of each site
for (s in 1:length(site_list)) { 
  # identify the input file-name based on the site_id
  siteFileNum    <- grep(site_list[s], sitesetFiles)
  inSiteForcFile <- sitesetFiles[siteFileNum]
  print(site_list[s])
  print(inSiteForcFile)
  fid <- nc_open(inSiteForcFile)
    lon_orig[s]<-ncvar_get(fid,"longitude")
    lat_orig[s]<-ncvar_get(fid,"latitude")
  nc_close(fid)
  # in case longitude needs to be converted back to the range [-180,180]
  if (lon_orig[s] > 180){ 
    lon_orig[s] <- lon_orig[s] - 360
    print(paste0("longitude converted to the range [-180,180], is now: ", lon_orig[s]))
  }
}
print("Done getting lon/lat from the metforcing.nc of each site")

#########################################################################################
## get supplementary data from Soil Grids Database
## new input path
print("")
print("getting supplementary data from Soil Grids Database")
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/"

# note:
# sand silt awc are read & used from the above site_list.dat file

ncid  <- nc_open(paste0(ipath,"SNDPPT.soilgrid.3600.1800.7.nc"))
lon_sg<-ncvar_get(ncid,'longitude')
lat_sg<-ncvar_get(ncid,'latitude')
sand_sg<-ncvar_get(ncid,'SNDPPT')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'SLTPPT.soilgrid.3600.1800.7.nc',sep=''))
silt_sg<-ncvar_get(ncid,'SLTPPT')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'BLDFIE.soilgrid.3600.1800.7.nc',sep=''))
bulk_density_sg<-ncvar_get(ncid,'BLDFIE')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'PHIHOX.soilgrid.3600.1800.7.nc',sep=''))
ph_sg<-ncvar_get(ncid,'PHIHOX')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'AWCtS.soilgrid.3600.1800.7.nc',sep=''))
fc_vol_sg<-ncvar_get(ncid,'AWCtS')
nc_close(ncid)

ncid<-nc_open(paste(ipath,'WWP.soilgrid.3600.1800.7.nc',sep=''))
pwp_vol_sg<-ncvar_get(ncid,'WWP')
nc_close(ncid)

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
ph           <- numeric(length(site_list))
bulk_density <- numeric(length(site_list))
sand         <- numeric(length(site_list))
silt         <- numeric(length(site_list))
awc          <- numeric(length(site_list))
awc_sl       <- array(NA,c(length(site_list),7))

## loop over sites: get soil data
for (s in 1:length(site_list)) {
  lonidx          <- which.min(abs(lon_sg-lon_orig[s]))
  latidx          <- which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  ph[s]           <- mean(ph_sg[lonidx,latidx,])/10
  depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]

  sand[s]     <- weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  silt[s]     <- weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  # clay is calculated below as: 1 - sand - silt
  awc[s]      <- weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
  awc_sl[s,]  <- (fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
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

ncid             <- nc_open(paste0(inPhosPath,"Labile_Inorganic_P.nc"))  # same lon/lat as above, reading it anyway to improve modularity & reduce dependencies 
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

  # if no value was found in this neares soil map location, extend the search
  if( any(is.na(c(soilP_labile_site[s],soilP_slow_site[s],soilP_occlud_site[s],soilP_primary_site[s]))) ){
    minDegree <- 0.8
    lonidx    <- which(abs(lon_sg-lon_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    latidx    <- which(abs(lat_sg-lat_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    # if at least one site could be found
    if (length(lonidx) > 0 && length(latidx) > 0){
      soilP_depth_site[s]   <- 0.5
      soilP_labile_site[s]  <- round(mean(soilP_labile_sg[lonidx,latidx],na.rm=T),4)
      soilP_slow_site[s]    <- round(mean(soilP_slow_sg[lonidx,latidx],na.rm=T),4)
      soilP_occlud_site[s]  <- round(mean(soilP_occlud_sg[lonidx,latidx],na.rm=T),4)
      soilP_primary_site[s] <- round(mean(soilP_primary_sg[lonidx,latidx],na.rm=T),4)
    }
  }

  # debug output
  #print(paste0("site: ",s," ",soilP_labile_site[s]," ",soilP_slow_site[s]," ",soilP_occlud_site[s]," ",soilP_primary_site[s]))
}

## replace potential NA values by the median of all values
soilP_labile_site[is.na(soilP_labile_site)]   <- round(median(soilP_labile_site, na.rm=T),4)
soilP_slow_site[is.na(soilP_slow_site)]       <- round(median(soilP_slow_site, na.rm=T),4)
soilP_occlud_site[is.na(soilP_occlud_site)]   <- round(median(soilP_occlud_site, na.rm=T),4)
soilP_primary_site[is.na(soilP_primary_site)] <- round(median(soilP_primary_site, na.rm=T),4)
soilP_depth_site[is.na(soilP_depth_site)]     <- round(median(soilP_depth_site, na.rm=T),4)


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
  if(!is.na(nwrbCat_sg[lonidx,latidx])){
    qmaxorg_site[s] <- qmaxorgTable$qmax_org_value[ as.integer(nwrbCat_sg[lonidx,latidx]) ]
  }
  # if no value was found in this neares soil map location, extend the search
  else{
    minDegree       <- 0.2
    lonidx          <- which(abs(lon_sg-lon_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    latidx          <- which(abs(lat_sg-lat_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    qmaxorg_site[s] <- qmaxorgTable$qmax_org_value[ as.integer(round(median(nwrbCat_sg[lonidx,latidx],na.rm=T),0)) ]
  }
  
  # round to 4 decimal places
  qmaxorg_site[s] <- round(qmaxorg_site[s],4)
  
  # debug output
  #print(paste0("site: ",s," ",qmaxorg_site[s]))
}

#########################################################################################
## output
## create empty output array and paste data into it
siteParam <- 28
siteinfo  <- array(NA,c(length(site_list),siteParam))

siteinfo[,1]  <- site_list
siteinfo[,2]  <- lon_orig
siteinfo[,3]  <- lat_orig
siteinfo[,4]  <- pft
siteinfo[,5]  <- allsitestartyear
siteinfo[,6]  <- allsiteendyear
siteinfo[,7]  <- lon_orig
siteinfo[,8]  <- clay
siteinfo[,9]  <- silt
siteinfo[,10] <- sand
siteinfo[,11] <- awc
siteinfo[,12] <- bulk_density
siteinfo[,13] <- ph
siteinfo[,14] <- taxousda
siteinfo[,15] <- taxnwrb
siteinfo[,16] <- glim_class_site
siteinfo[,17] <- lai
siteinfo[,18] <- nleaf
siteinfo[,19] <- sla
siteinfo[,20] <- height
siteinfo[,21] <- plantyear
siteinfo[,22] <- soilP_depth_site
siteinfo[,23] <- soilP_labile_site
siteinfo[,24] <- soilP_slow_site
siteinfo[,25] <- soilP_occlud_site
siteinfo[,26] <- soilP_primary_site
siteinfo[,27] <- qmaxorg_site
siteinfo[,28] <- allsiteorigin

## output
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp','origin')

# 1991-2018  (the length of the available forcing data varies across fluxnet sites, but lies within 1991-2018)
write.table(siteinfo,file=paste0("plumber2_all_sites_1991-2018_list",Sys.Date(),".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1901
# 1901-... end year same as above
write.table(siteinfo,file=paste0("plumber2_all_sites_1901-2018_list",Sys.Date(),".dat"),quote=F,row.names=F,col.names=T)

print("")
print(paste0("Script has finished. The two new plumber2_all_sites_*.dat files with ",siteParam," parameters were created."))
print(" Please remove the DATE in the file name and copy the two files to /Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET3/ ")


## end ##
