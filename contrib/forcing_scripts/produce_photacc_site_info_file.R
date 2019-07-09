##
## generate site_list.dat files for the site set PHOTACC (phot acclimation), which are cruncep type of sites ticket #557 for details
##

## library
library(ncdf4)

#########################################################################################
## run this script by:
## cd contrib/forcing_scripts/
## Rscript --vanilla ./produce_photacc_site_info_file.R


#########################################################################################
## output 
## the output files will be written in this folder


#########################################################################################
## basic site lists
metsites      <- read.table("generate_photacc_siteinfo_noheader.dat", stringsAsFactors=F, header=F)
metsite_list  <- metsites[,1]


#########################################################################################
## available site characteristics
dat           <- metsites
colnames(dat) <- c("Site-ID", "lon", "lat", "PFT")

# number of sites
siteNum <- nrow(metsites)

## PFT information per site
pft <- dat$PFT

## more site characteristics
siteID_list     <- dat$'Site-ID'
lon_orig        <- dat$lon
lat_orig        <- dat$lat
plantyear       <- c(rep(-9999,siteNum))

## feedback
print(paste0(" Start creating the site_list info file for ",siteNum," sites of the PHOTACC siteset."))


#########################################################################################
## get supplementary data from Soil Grids Database
## input path
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/"

ncid<-nc_open(paste(ipath,'SNDPPT.soilgrid.3600.1800.7.nc',sep=''))
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
taxnwrb      <- numeric(siteNum)
taxousda     <- numeric(siteNum)
depth        <- numeric(siteNum)
ph           <- numeric(siteNum)
bulk_density <- numeric(siteNum)
sand         <- numeric(siteNum)
silt         <- numeric(siteNum)
awc_sl       <- array(NA,c(siteNum,7))
awc          <- numeric(siteNum)

for (s in 1:siteNum) {
  lonidx          <- which.min(abs(lon_sg-lon_orig[s]))
  latidx          <- which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  ph[s]           <- mean(ph_sg[lonidx,latidx,])/10
  depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]

  sand[s]     <- weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  silt[s]     <- weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  awc_sl[s,]  <- (fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  awc[s]      <- weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}


#########################################################################################
## read lithology map
## new input file
ncid <- nc_open("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc")
lon_glim        <- ncvar_get(ncid,'lon')
lat_glim        <- ncvar_get(ncid,'lat')
glim_class      <- ncvar_get(ncid,'GLim')
glim_class_site <- numeric(siteNum)
nc_close(ncid)

## loop over sites: get soil data
for (s in 1:siteNum) {
  lonidx             <- which.min(abs(lon_glim-lon_orig[s]))
  latidx             <- which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s] <- glim_class[lonidx,latidx]
}


#########################################################################################
## replace NAN with default values and run some additional calculations
sand[is.na(sand)] <- 0.4
silt[is.na(silt)] <- 0.4
clay              <- 1 - sand - silt
ph[is.na(ph)]     <- 6
awc[is.na(awc)]   <- 200
plantyear[is.na(plantyear)]             <- 1500
bulk_density[is.na(bulk_density)]       <- 1500
taxousda[is.na(taxousda)]               <- 30
taxnwrb[is.na(taxnwrb)]                 <- 27
glim_class_site[is.na(glim_class_site)] <- 2



#########################################################################################
## read site-specific data on Phosphorous
## multiple new input files
inPhosPath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/Phosphorous/v2014_06/Data/"

# soilP_depth_sg is 0.5 [m], same value for all sites

ncid             <- nc_open(paste0(inPhosPath,"Labile_Inorganic_P.nc"))   # same lon/lat as above, reading it anyway to improve modularity & reduce dependencies 
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
soilP_depth_site   <- numeric(siteNum)
soilP_labile_site  <- numeric(siteNum)
soilP_slow_site    <- numeric(siteNum)
soilP_occlud_site  <- numeric(siteNum)
soilP_primary_site <- numeric(siteNum)

## loop over sites: get soil P data from nearest location
for (s in 1:siteNum) {
  # get neares lon/lat in soil map
  lonidx                <- which.min(abs(lon_sg-lon_orig[s]))
  latidx                <- which.min(abs(lat_sg-lat_orig[s]))
  
  soilP_depth_site[s]   <- 0.5
  soilP_labile_site[s]  <- round(soilP_labile_sg[lonidx,latidx],4)
  soilP_slow_site[s]    <- round(soilP_slow_sg[lonidx,latidx],4)
  soilP_occlud_site[s]  <- round(soilP_occlud_sg[lonidx,latidx],4)
  soilP_primary_site[s] <- round(soilP_primary_sg[lonidx,latidx],4)
  

  # all sites get a value by the above "nearest neighbour" approach (2018-05-16), so the "extended search" is not needed 
  #
  # # if no value was found in this neares soil map location, extend the search
  # if( any(is.na(c(soilP_labile_site[s],soilP_slow_site[s],soilP_occlud_site[s],soilP_primary_site[s]))) ){
  #   minDegree <- 0.8
  #   lonidx    <- which(abs(lon_sg-lon_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
  #   latidx    <- which(abs(lat_sg-lat_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
  #   # if at least one site could be found
  #   if (length(lonidx) > 0 && length(latidx) > 0){
  #     soilP_depth_site[s]   <- 0.5
  #     soilP_labile_site[s]  <- round(mean(soilP_labile_sg[lonidx,latidx],na.rm=T),4)
  #     soilP_slow_site[s]    <- round(mean(soilP_slow_sg[lonidx,latidx],na.rm=T),4)
  #     soilP_occlud_site[s]  <- round(mean(soilP_occlud_sg[lonidx,latidx],na.rm=T),4)
  #     soilP_primary_site[s] <- round(mean(soilP_primary_sg[lonidx,latidx],na.rm=T),4)
  #   }
  # }
  
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
qmaxorg_site <- numeric(siteNum)

## loop over sites: get NWRB cat data from nearest location & get according Qmax_org
for (s in 1:siteNum) {
  # get neares lon/lat in soil map
  lonidx <- which.min(abs(lon_sg-lon_orig[s]))
  latidx <- which.min(abs(lat_sg-lat_orig[s]))
  if(!is.na(nwrbCat_sg[lonidx,latidx])){
    qmaxorg_site[s] <- qmaxorgTable$qmax_org_value[ as.integer(nwrbCat_sg[lonidx,latidx]) ]
  }
  # if no value was found in this neares soil map location, extend the search
  else{
    minDegree       <- 0.4
    lonidx          <- which(abs(lon_sg-lon_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    latidx          <- which(abs(lat_sg-lat_orig[s]) < minDegree) # all points within a distance of less than minDegree degree
    qmaxorg_site[s] <- qmaxorgTable$qmax_org_value[ as.integer(round(median(nwrbCat_sg[lonidx,latidx],na.rm=T),0)) ]
  }
  
  # round to 4 decimal places
  qmaxorg_site[s] <- round(qmaxorg_site[s],4)
  
  # debug output
  #print(paste0("site: ",s," ",qmaxorg_site[s]))
}

## if no value for a site is found, setting it to a default value (using simply the median of all sites)
qmaxorg_site[is.na(qmaxorg_site)] <- round(median(qmaxorg_site, na.rm=T),4)


#########################################################################################
## start creating the final output file: photacc_sites_1901-2015
##  

# output file
outfileName   <- "photacc_sites_1901-2015"
siteParamNum  <- 27
siteinfo      <- array(NA,c(siteNum, siteParamNum))
startYear     <- 1901
lastYear      <- 2015

# feedback
print(" All data are gathered.")
print(paste0(" Start writing the file ",outfileName,"DATE.dat. Please remove DATE before using this file."))

# loop over all photacc sites
for (i in 1:siteNum) {
    siteinfo[i,1]  <- siteID_list[i]
    siteinfo[i,2]  <- lon_orig[i]
    siteinfo[i,3]  <- lat_orig[i]
    siteinfo[i,4]  <- pft[i]
    siteinfo[i,5]  <- startYear
    siteinfo[i,6]  <- lastYear
    siteinfo[i,7]  <- lon_orig[i]
    siteinfo[i,8]  <- clay[i]
    siteinfo[i,9]  <- silt[i]
    siteinfo[i,10] <- sand[i]
    siteinfo[i,11] <- awc[i]
    siteinfo[i,12] <- bulk_density[i]
    siteinfo[i,13] <- ph[i]
    siteinfo[i,14] <- taxousda[i]
    siteinfo[i,15] <- taxnwrb[i]
    siteinfo[i,16] <- glim_class_site[i]
    siteinfo[i,17] <- -9999   #LAI
    siteinfo[i,18] <- -9999   #Nleaf
    siteinfo[i,19] <- -9999   #SLA
    siteinfo[i,20] <- -9999   #Height
    siteinfo[i,21] <- plantyear[i]
    siteinfo[i,22] <- soilP_depth_site[i]
    siteinfo[i,23] <- soilP_labile_site[i]
    siteinfo[i,24] <- soilP_slow_site[i]
    siteinfo[i,25] <- soilP_occlud_site[i]
    siteinfo[i,26] <- soilP_primary_site[i]
    siteinfo[i,27] <- qmaxorg_site[i]
}

## output
# header (must match siteParamNum)
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')
# 1901-2015 data
write.table(siteinfo, file=paste0(outfileName, Sys.Date(), ".dat"), quote=F, row.names=F, col.names=T)


# feedback
print(" Done!")


## fin ##
