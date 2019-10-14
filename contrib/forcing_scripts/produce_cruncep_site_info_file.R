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
## generate site_list.dat files for cruncep site sets for QUINCY (the spp_gfdb may be used for SPP_1685 site-set too): 
##    all_gfdb  cue_gfdb  spp_gfdb
##

## library
library(ncdf4)

#########################################################################################
## run this script by:
## cd contrib/forcing_scripts/
## Rscript --vanilla ./produce_cruncep_site_info_file.R


#########################################################################################
## output 
## the output files will be written in this folder


#########################################################################################
## basic input path
baseInputPath <- "/Net/Groups/BGI/people/szaehle/Projects/QUINCY/GenerateCRUForcing_2017/"


#########################################################################################
## basic site lists
## (orininally from: baseInputPath,'siteinfo_noheader.dat' & 'addsites2019_1.dat', moved into repo for reproducability)
# all & cue cruncep
metsites      <- read.table("generate_cruncep_siteinfo_noheader.dat", stringsAsFactors=F, header=F)
metsite_list  <- metsites[,1]
# spp cruncep
metsites_spp  <- read.table("generate_cruncep_spp_siteinfo_noheader.dat", stringsAsFactors=F, header=F)
spp_site_list <- metsites_spp[,1]


#########################################################################################
## available site characteristics (incudes multiple "sub-sites" per site, e.g. AND_01 AND_02 AND_03 ..)
dat<-read.csv(paste0(baseInputPath,'info_forest.csv'),sep=';',stringsAsFactors=F)

# number of sites in the the info_forest.csv; last site in info_forest.csv is ZAF with some more rows after this site with no information
forestSiteNum <- which(dat[,1] == "ZAF") # should be 580 (after which empty data-rows follow in this file)

## PFT information per site
# create pft vector from information of pft1 pft2 pft3
pft1 <- numeric(nrow(dat))
pft2 <- numeric(nrow(dat))
pft3 <- numeric(nrow(dat))
# pft1
pft1[dat[,3]=='Boreal Humid']='B'
pft1[dat[,3]=='Boreal Semi-arid']='B'
pft1[dat[,3]=='Mediterranean Warm']='Te'
pft1[dat[,3]=='Temperate Humid']='Te'
pft1[dat[,3]=='Temperate Semi-arid']='Te'
pft1[dat[,3]=='Tropical Humid']='Tr'
pft1[dat[,3]=='Tropical Semi-arid']='Tr'
# pft2 
pft2[dat[,4]=='Broadleaved']='B'
pft2[dat[,4]=='Mixed']='M'
pft2[dat[,4]=='Needle-leaved']='N'
# pft3
pft3[dat[,5]=='Deciduous']='S'
pft3[dat[,5]=='Mixed']='M'
pft3[dat[,5]=='Evergreen']='E'
pft3[pft1=='Tr'&pft3=='S']<-'R'
# change pft1
pft1[pft1=='Te'&pft2=='N'&pft3=='S']<-'B'
pft1[pft1=='Tr'&pft2=='N'&pft3=='E']<-'Te'
# combine pft1,2,3 into the pft vector
pft <- paste0(pft1,pft2,pft3)

## more site characteristics
forestsite_list <- dat[,1]
lon_orig        <- dat$Longitude
lat_orig        <- dat$Latitude
sand            <- dat$Sand /100
silt            <- dat$Silt /100
awc             <- dat$Profile_avail_WATER
age             <- dat$Stand.age
nstat           <- dat$NutrientAvailability
nstat[is.na(nstat)]<-0


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
taxnwrb      <- numeric(forestSiteNum)
taxousda     <- numeric(forestSiteNum)
depth        <- numeric(forestSiteNum)
ph           <- numeric(forestSiteNum)
bulk_density <- numeric(forestSiteNum)
awc_sl       <- array(NA,c(forestSiteNum,7))

## loop over all sites in info_forest.csv and get soil data from above read netcdf files
print(paste0("Creating site info file for selected sites of the ",forestSiteNum," cruncep sites from the info_forest.csv."))
print(" ")
for (s in 1:forestSiteNum) {
  lonidx          <- which.min(abs(lon_sg-lon_orig[s]))
  latidx          <- which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  ph[s]           <- mean(ph_sg[lonidx,latidx,])/10
  depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]

  if(is.na(sand[s])) sand[s] <- weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  if(is.na(silt[s])) silt[s] <- weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  awc_sl[s,]                 <-(fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  if(is.na(awc[s]))  awc[s]  <- weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}


#########################################################################################
## read lithology map
## new input file
ncid <- nc_open("/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc")
lon_glim        <- ncvar_get(ncid,'lon')
lat_glim        <- ncvar_get(ncid,'lat')
glim_class      <- ncvar_get(ncid,'GLim')
glim_class_site <- numeric(forestSiteNum)
nc_close(ncid)

## loop over sites: get soil data
for (s in 1:forestSiteNum) {
  lonidx             <- which.min(abs(lon_glim-lon_orig[s]))
  latidx             <- which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s] <- glim_class[lonidx,latidx]
}


#########################################################################################
## replace NAN with default values
sand[is.na(sand)] <- 0.4
silt[is.na(silt)] <- 0.4
clay              <- 1 - sand - silt
ph[is.na(ph)]     <- 6
awc[is.na(awc)]   <- 200
age[is.na(age)]   <- 500
plantyear         <- 2000 - age
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
soilP_depth_site   <- numeric(forestSiteNum)
soilP_labile_site  <- numeric(forestSiteNum)
soilP_slow_site    <- numeric(forestSiteNum)
soilP_occlud_site  <- numeric(forestSiteNum)
soilP_primary_site <- numeric(forestSiteNum)

## loop over sites: get soil P data from nearest location
for (s in 1:forestSiteNum) {
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
qmaxorg_site <- numeric(forestSiteNum)

## loop over sites: get NWRB cat data from nearest location & get according Qmax_org
for (s in 1:forestSiteNum) {
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

## no value for site number 434 "SUE", so setting it to a default value (using simply the median of all sites)
# to be save, set all sites with NA to a default value
qmaxorg_site[is.na(qmaxorg_site)] <- round(median(qmaxorg_site, na.rm=T),4)



#########################################################################################
## get all gfdb sites that match specific criteria and create a vector with consistent siteID naming scheme
##    selected sites:    gfdb_all_selected_sitelist   (subset of forestsite_list)
##    consistent siteID: gfdb_sitename_new
## use available site information to constrain sites: GPP, NPP and sections of it  (even if not matching)
## new input file (in addition to the above read file: info_forest.csv)
gfdb_all <- read.csv(paste0(baseInputPath,'info_forest_DATA.csv'),sep=';',stringsAsFactors=F)

## get True/False vector for all forestSiteNum (580) sites
# where True = sites that have values for TNPP, GPP and ANPP but are not the individually named sites such as "MAN"
# (sort out sites with no data for TNPP, GPP or ANPP)
all_idx <- !(is.na(gfdb_all$TNPP)&is.na(gfdb_all$GPP)&is.na(gfdb_all$ANPP))&
           gfdb_all[,1]!='MAN'&gfdb_all[,1]!='PUU'&gfdb_all[,1]!='SKO'&gfdb_all[,1]!='SUE'&
           gfdb_all[,1]!='YUS_01'&gfdb_all[,1]!='YUS_02'&gfdb_all[,1]!='CHAM_01'&
           gfdb_all[,1]!='CHAM_02'&gfdb_all[,1]!='CHAM_03'

# all sites marked as True
dummy_selected          <- gfdb_all[all_idx,]
dummy_selected_sitelist <- dummy_selected[,1]
not_mixed               <- logical(length(dummy_selected_sitelist))  # sets all to False
i <- 1

## loop over all "True" sites: set the vector not_mixed to True for sites that have no mixed PFTs, i.e., just one PFT
# (sort out sites with more than one PFT)
for (gfdb_site in dummy_selected_sitelist) {
  # get a True/False vector with with length of forestsite_list that is True at the position of gfdb_site ID
  idx <- gfdb_site==forestsite_list    # siteID from info_forest.csv
  if(pft2[idx]!='M'&pft3[idx]!='M'){
    not_mixed[i]<-TRUE}
  else{
    not_mixed[i]<-FALSE}
  i<-i+1
}

## vector with the selected sites
gfdb_all_selected          <- dummy_selected[not_mixed,]
gfdb_all_selected_sitelist <- gfdb_all_selected[,1]


## correct name convention by matching GFDB data with metsite database
oldsite            <- "" 
newsite            <- "" 
metsite_list_found <- c()
gfdb_sitename_new  <- numeric(length(gfdb_all_selected_sitelist))
count <- 1
i     <- 1
# loop over the selected sites: correct site name
for (gfdb_site in gfdb_all_selected_sitelist){
 dummy <- gfdb_site
 # renamed site list
 if(dummy=='BILY') dummy<-'BLY'
 if(dummy=='CASTEL') dummy<-'CAP'
 if(dummy=='CHAK') dummy<-'CHX'
 if(strtrim(dummy,4)=='CHAM') dummy<-'CHM'
 if(dummy=='FRAN') dummy<-'FRN'
 if(strtrim(dummy,4)=='HARD') dummy<-'HRD'
 if(dummy=='HEST') dummy<-'HST'
 if(dummy=='KOIN') dummy<-'KIN'
 if(dummy=='KOIW') dummy<-'KIW'
 if(dummy=='LAUPA') dummy<-'LPA'
 if(dummy=='LOOB') dummy<-'LOB'
 if(dummy=='PRI_SOAS') dummy<-'PRS'
 if(dummy=='PRI_SOBS') dummy<-'PRB'
 if(dummy=='PRI_SOJP') dummy<-'PRP'
 if(dummy=='PRI_SYJP') dummy<-'PRP'
 if(dummy=='THO_S') dummy<-'THS'
 if(dummy=='TOMA') dummy<-'TMA'
 if(dummy=='TURA') dummy<-'TRA'
 if(dummy=="THO_W13") {oldsite="THO";count=6} 
 if(dummy=='USA_01') dummy<-'US1'
 if(dummy=='USA_02') dummy<-'US2'
 if(dummy=='USA_03') dummy<-'US2'
 if(dummy=='USA_04') dummy<-'US2'
 if(dummy=='USA_05') dummy<-'US2'
 if(dummy=='WAL_01') dummy<-'WL0'
 if(dummy=='WAL_02') dummy<-'WL0'
 if(dummy=='WAL_03') dummy<-'WL0'
 if(dummy=='WAL_04') dummy<-'WL0'
 if(dummy=='WC_00') dummy<-'WCA'

 # check if this site is available in metsite_list | met forcing data file read above: siteinfo_noheader.dat
 idx<-metsite_list==strtrim(dummy,3)
 if(any(idx)) {
    newsite <- metsite_list[idx]
    print(paste0("gfdb selected site: ",dummy," matches with the met forcing site: ",newsite))
    
    # compare newsite name with oldsite name to account for multiple "sub-sites" per site, e.g. AND_01 AND_02 AND_03 ..
    # increase count only in case the site names match
    if(newsite == oldsite) {count <- count+1}
    else {count=1}
    if(count<10)  {gfdb_sitename_new[i] <- paste0(metsite_list[idx],"_0",count)}
    if(count>=10) {gfdb_sitename_new[i] <- paste0(metsite_list[idx],"_",count)}

    oldsite <- newsite
    i<-i+1
 }
 else {print(paste0("|| the gfdb selected site: ",dummy," does NOT match any met forcing site"))}
}


#########################################################################################
## start creating the final output file: all_gfdb_sites* 1901-2015 and 1984-2013
##  
## using the objects:
##  forestsite_list             -> list of siteID from info_forest.csv (sites specific data such as PFT lon/lat soil-composition)
##                                 -> this list is used to generate all the data vectors, so use this for "siteinfo" creation 
##  gfdb_all_selected_sitelist  -> the selected gfdb sites (matching data availability criteria)
##                                 -> this is the list of valid site names
##  gfdb_cue_sitename_new       -> gfdb site names (ID) in output file (matching the metsite database)

# output file
siteParam <- 27
siteinfo  <- array(NA,c(length(gfdb_all_selected_sitelist),siteParam))
# additional output file with original gdfb data (info_forest_DATA.csv) for selected sites only
sitedata  <- array(NA,c(length(gfdb_all_selected_sitelist),length(gfdb_all_selected[1,])+1))
# loop over gfdb selected sites
i<-1
for (gfdb_site in gfdb_all_selected_sitelist) {
  idx  <- gfdb_site==forestsite_list             # all sites | siteID from info_forest.csv
  idx2 <- gfdb_site==gfdb_all_selected_sitelist  # selected sites with | identical siteID as in forestsite_list (info_forest.csv)
    siteinfo[i,1]<-gfdb_sitename_new[idx2]       # renamed siteID of selected sites (all siteID are: XXX_NN [X name; N number])
    siteinfo[i,2]<-lon_orig[idx]
    siteinfo[i,3]<-lat_orig[idx]
    siteinfo[i,4]<-pft[idx]
    siteinfo[i,5]<-1901
    siteinfo[i,6]<-2015
    siteinfo[i,7]<-lon_orig[idx]
    siteinfo[i,8]<-clay[idx]
    siteinfo[i,9]<-silt[idx]
    siteinfo[i,10]<-sand[idx]
    siteinfo[i,11]<-awc[idx]
    siteinfo[i,12]<-bulk_density[idx]
    siteinfo[i,13]<-ph[idx]
    siteinfo[i,14]<-taxousda[idx]
    siteinfo[i,15]<-taxnwrb[idx]
    siteinfo[i,16]<-glim_class_site[idx]
    siteinfo[i,17]<-(-9999) #LAI
    siteinfo[i,18]<-(-9999) #Nleaf
    siteinfo[i,19]<-(-9999) #SLA
    siteinfo[i,20]<-(-9999) #Height
    siteinfo[i,21]<-plantyear[idx]
    siteinfo[i,22]<-soilP_depth_site[idx]
    siteinfo[i,23]<-soilP_labile_site[idx]
    siteinfo[i,24]<-soilP_slow_site[idx]
    siteinfo[i,25]<-soilP_occlud_site[idx]
    siteinfo[i,26]<-soilP_primary_site[idx]
    siteinfo[i,27]<-qmaxorg_site[idx]
    
    
    # extra output
    sitedata[i,1]<-gfdb_sitename_new[idx2]
    for (j in 2:11) sitedata[i,j]<-gfdb_all_selected[idx2,j]
    sitedata[i,12]<-gfdb_all_selected[idx2,1]
    i<-i+1
}

## output
# time stamp
date <- Sys.Date()
# header
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')
# 1901-2015
write.table(siteinfo,file=paste0("crun_all_gfdb_sites_1901-2015_list",date,".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
# 1984-2013
write.table(siteinfo,file=paste0("crun_all_gfdb_sites_1984-2013_list",date,".dat"),quote=F,row.names=F,col.names=T)

## extra output with specific date
colnames(sitedata)<-c(colnames(gfdb_all_selected),'Original.Site.Code')
write.table(sitedata,file=paste0("global_forest_database_all_sites_",date,".csv"),quote=F,row.names=F,col.names=T,sep=",")



#########################################################################################
## preparing & creating the output file: cue_gfdb_sites* 1901-2015 and 1984-2013
##  
## repeating some site-selection steps as above but limited to cue sites of the gfdb sites
## no new data (for all_site_list.dat) is generated here
##

## new input file with available site information (CUE and compontents)
#  additionally to the above read file: info_forest.csv
#  but different from the above: gfdb_all <- read.csv('info_forest_DATA.csv'...) input file
gfdb_cue <- read.csv(paste0(baseInputPath,'global_forest_database_21092016.csv'),sep=',',stringsAsFactors=F)

all_idx                 <-gfdb_cue[,1]!='XIS_01'&gfdb_cue[,1]!='XIS_02'
dummy_selected          <-gfdb_cue[all_idx,]
dummy_selected_sitelist <-dummy_selected[,1]
not_mixed               <-logical(length(dummy_selected_sitelist))
i<-1

## loop over all intially selected sites 
#  and sort out sites with more than one PFT (mixed PFTs)
for (gfdb_site in dummy_selected_sitelist) {
  idx<-gfdb_site==forestsite_list    # list of siteID from info_forest.csv
  if(gfdb_site=='CAL')     idx<-'REN'==forestsite_list
  if(gfdb_site=='NAU_M')   idx<-'NAU'==forestsite_list
  if(gfdb_site=='NAU_U')   idx<-'NAU'==forestsite_list
  if(gfdb_site=='POPF')    idx<-'POP_01'==forestsite_list
  if(gfdb_site=='THO_D15') idx<-'THO_D12'==forestsite_list
  if(gfdb_site=='THO_D71') idx<-'THO_71'==forestsite_list
  if(pft2[idx]!='M'&pft3[idx]!='M') {not_mixed[i]<-TRUE}
  else {not_mixed[i]<-FALSE}
  i<-i+1
}

## vector with the selected sites
gfdb_cue<-dummy_selected[not_mixed,]


## correct name convention by matching GFDB data with metsite database
oldsite               <- "" 
newsite               <- "" 
gfdb_cue_sitelist     <- gfdb_cue[,1]
gfdb_cue_sitename_new <- numeric(length(gfdb_cue_sitelist))
count <- 1
i     <- 1
# loop over the selected sites: correct site name
for (gfdb_site in gfdb_cue_sitelist){
 dummy <- gfdb_site
 # renamed site list
 if(dummy=='BILY') dummy<-'BLY'
 if(dummy=='CAL') dummy<-'REN' #only for met purposes!!
 if(dummy=='CASTEL') dummy<-'CAP'
 if(dummy=='CHAK') dummy<-'CHX'
 if(strtrim(dummy,4)=='CHAM') dummy<-'CHM'
 if(dummy=='FRAN') dummy<-'FRN'
 if(strtrim(dummy,4)=='HARD') dummy<-'HRD'
 if(dummy=='HEST') dummy<-'HST'
 if(dummy=='KOIN') dummy<-'KIN'
 if(dummy=='KOIW') dummy<-'KIW'
 if(dummy=='LAUPA') dummy<-'LPA'
 if(dummy=='LOOB') dummy<-'LOB'
 if(dummy=='PRI_SOAS') dummy<-'PRS'
 if(dummy=='PRI_SOBS') dummy<-'PRB'
 if(dummy=='PRI_SOJP') dummy<-'PRP'
 if(dummy=='PRI_SYJP') dummy<-'PRP'
 if(dummy=='THO_S') dummy<-'THS'
 if(dummy=='TOMA') dummy<-'TMA'
 if(dummy=='TURA') dummy<-'TRA'
 if(dummy=="THO_W13") {oldsite="THO";count=7} 
 if(dummy=='USA_01') dummy<-'US1'
 if(dummy=='USA_02') dummy<-'US2'
 if(dummy=='USA_03') dummy<-'US2'
 if(dummy=='USA_04') dummy<-'US2'
 if(dummy=='USA_05') dummy<-'US2'
 if(dummy=='WAL_01') dummy<-'WL0'
 if(dummy=='WAL_02') dummy<-'WL0'
 if(dummy=='WAL_03') dummy<-'WL0'
 if(dummy=='WAL_04') dummy<-'WL0'
 if(dummy=='WC_00') dummy<-'WCA'

 # check if this site is available in metsite_list | met forcing data file read above: siteinfo_noheader.dat
 idx <- metsite_list==strtrim(dummy,3)
 if(any(idx)) {
    newsite <- metsite_list[idx]
    print(paste0("gfdb cue selected site: ",dummy," matches with the met forcing site: ",newsite))

    # compare newsite name with oldsite name to account for multiple "sub-sites" per site, e.g. AND_01 AND_02 AND_03 ..
    # increase count only in case the site names match
    if(newsite==oldsite) {count=count+1}
    else {count=1}
    if(count<10) {gfdb_cue_sitename_new[i]<-paste(metsite_list[idx],'_0',count,sep='')}
    if(count>=10) {gfdb_cue_sitename_new[i]<-paste(metsite_list[idx],'_',count,sep='')}
    
    oldsite <- newsite
    i<-i+1
  } 
 else {print(paste0("|| the gfdb cue selected site: ",dummy," does NOT match any met forcing site"))}
}


#########################################################################################
## prepare the output file: cue_gfdb_sites* 1901-2015 and 1984-2013
##  

# output file
siteParam <- 27
siteinfo  <- array(NA,c(length(gfdb_cue_sitelist),siteParam))
sitedata  <- array(NA,c(length(gfdb_cue_sitelist),length(gfdb_cue[1,])+1))
# loop over gfdb cue selected sites
i<-1
for (gfdb_site in gfdb_cue_sitelist) {
  idx <- gfdb_site==forestsite_list    # list of siteID from info_forest.csv
  if(gfdb_site=='CAL')     idx<-'REN'==forestsite_list
  if(gfdb_site=='NAU_M')   idx<-'NAU'==forestsite_list
  if(gfdb_site=='NAU_U')   idx<-'NAU'==forestsite_list
  if(gfdb_site=='POPF')    idx<-'POP_01'==forestsite_list
  if(gfdb_site=='THO_D15') idx<-'THO_D12'==forestsite_list
  if(gfdb_site=='THO_D71') idx<-'THO_71'==forestsite_list
  idx2<-gfdb_site==gfdb_cue_sitelist
  
    siteinfo[i,1]<-gfdb_cue_sitename_new[idx2]
    siteinfo[i,2]<-lon_orig[idx]
    siteinfo[i,3]<-lat_orig[idx]
    siteinfo[i,4]<-pft[idx]
    siteinfo[i,5]<-1901
    siteinfo[i,6]<-2015
    siteinfo[i,7]<-lon_orig[idx]
    siteinfo[i,8]<-clay[idx]
    siteinfo[i,9]<-silt[idx]
    siteinfo[i,10]<-sand[idx]
    siteinfo[i,11]<-awc[idx]
    siteinfo[i,12]<-bulk_density[idx]
    siteinfo[i,13]<-ph[idx]
    siteinfo[i,14]<-taxousda[idx]
    siteinfo[i,15]<-taxnwrb[idx]
    siteinfo[i,16]<-glim_class_site[idx]
    siteinfo[i,17]<-(-9999) #LAI
    siteinfo[i,18]<-(-9999) #Nleaf
    siteinfo[i,19]<-(-9999) #SLA
    siteinfo[i,20]<-(-9999) #Height
    siteinfo[i,21]<-plantyear[idx]
    siteinfo[i,22]<-soilP_depth_site[idx]
    siteinfo[i,23]<-soilP_labile_site[idx]
    siteinfo[i,24]<-soilP_slow_site[idx]
    siteinfo[i,25]<-soilP_occlud_site[idx]
    siteinfo[i,26]<-soilP_primary_site[idx]
    siteinfo[i,27]<-qmaxorg_site[idx]
    
    
    # extra output
    sitedata[i,1]<-gfdb_cue_sitename_new[idx2]
    for (j in 2:length(gfdb_cue[1,])) sitedata[i,j]<-gfdb_cue[idx2,j]
    sitedata[i,length(gfdb_cue[1,])+1]<-gfdb_cue[idx2,1]
    i<-i+1
}


## output
# time stamp
date <- Sys.Date()
# header
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')
# 1901-2015
write.table(siteinfo,file=paste0("crun_cue_gfdb_sites_1901-2015_list",date,".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
# 1984-2013
write.table(siteinfo,file=paste0("crun_cue_gfdb_sites_1984-2013_list",date,".dat"),quote=F,row.names=F,col.names=T)

## extra output with specific date
colnames(sitedata)<-c(colnames(gfdb_cue),'Original.Site.Code')
write.table(sitedata,file=paste0("global_forest_database_cue_sites_",date,".csv"),quote=F,row.names=F,col.names=T,sep=",")




#########################################################################################
## preparing & creating the output file: spp_gfdb_sites* 1901-2015 and 1984-2013
## 
## 
## new data is read from files / generated here
##
##

## get number of spp sites (see file input above)
sppSiteNum <- length(spp_site_list)

## get the lon/lat of all sites from the metsites_spp object (with list of siteID in metsites_spp), file siteinfo_noheader.dat
#  new values for lon_orig & lat_orig  (above they contain lon/lat of all sites in forestsite_list from info_forest.csv)
lon_orig <- numeric(sppSiteNum)
lat_orig <- numeric(sppSiteNum)
i<-1
for (spp_site in spp_site_list) {
    idx         <- spp_site==metsites_spp
    lon_orig[i] <- metsites_spp[idx,2]
    lat_orig[i] <- metsites_spp[idx,3]
  i<-i+1
}


## update the lon/lat values of the soilgrid used below for deriving soil data
#  the data objects, read from nc files, should still be available here (i.e. not removed in the all_cruncep section above)
#
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/"

ncid   <- nc_open(paste0(ipath,"SNDPPT.soilgrid.3600.1800.7.nc"))
lon_sg <- ncvar_get(ncid,'longitude')
lat_sg <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create new (empty) vectors/arrays for output
awc          <- numeric(sppSiteNum)
clay         <- numeric(sppSiteNum)
silt         <- numeric(sppSiteNum)
sand         <- numeric(sppSiteNum)
taxnwrb      <- numeric(sppSiteNum)
taxousda     <- numeric(sppSiteNum)
depth        <- numeric(sppSiteNum)
ph           <- numeric(sppSiteNum)
bulk_density <- numeric(sppSiteNum)
awc_sl       <- array(NA,c(sppSiteNum,7))

## loop over the sites and get soil data from above read netcdf files
for (s in 1:sppSiteNum) {
  lonidx <- which.min(abs(lon_sg-lon_orig[s]))
  latidx <- which.min(abs(lat_sg-lat_orig[s]))
  
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]
  awc_sl[s,]      <- (fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  awc[s]          <- weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
  
  #ph[s]           <- mean(ph_sg[lonidx,latidx,])/10
  #sand[s]         <- weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  #silt[s]         <- weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
}
## set remaining parameter from known data, i.e., Lin Yu provided these data (proportions of soil types and ph)
# site order: BBR LUE VES MTF CON
sand    <- c(22.59, 76.6,  39.0, 51.0, 41.0) / 100.0
silt    <- c(51.31, 17.13, 45.0, 31.0, 42.0) / 100.0
clay    <- 1.0 - sand - silt
ph      <- c(5.2, 4.3, 4.34, 4.46, 4.42)
pft_spp <- c("TeBS", "TeBS", "TeBS", "TeBS", "TeBS")

## read lithology map
#  this is the same input file as above, but repeating it here again for improved modularity
#
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/"

ncid            <- nc_open(paste0(ipath,"GLim.720.360.nc"))
lon_glim        <- ncvar_get(ncid,'lon')
lat_glim        <- ncvar_get(ncid,'lat')
glim_class      <- ncvar_get(ncid,'GLim')
glim_class_site <- numeric(sppSiteNum)
nc_close(ncid)

## loop over sites: get soil data
for (s in 1:sppSiteNum) {
  lonidx             <- which.min(abs(lon_glim-lon_orig[s]))
  latidx             <- which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s] <- glim_class[lonidx,latidx]
}


## read site-specific data on Phosphorous
#  only update lon_sg/lat_sg but use the data objects alread read above (all_cruncep section)
#
inPhosPath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/Phosphorous/v2014_06/Data/"

ncid             <- nc_open(paste0(inPhosPath,"Labile_Inorganic_P.nc"))
lon_sg           <- ncvar_get(ncid,'longitude')
lat_sg           <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create empty vectors for site-data
soilP_depth_site   <- numeric(sppSiteNum)
soilP_labile_site  <- numeric(sppSiteNum)
soilP_slow_site    <- numeric(sppSiteNum)
soilP_occlud_site  <- numeric(sppSiteNum)
soilP_primary_site <- numeric(sppSiteNum)

## loop over the two sites: get soil P data from nearest location
for (s in 1:sppSiteNum) {
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


## get site-specific data on qmax_org_fine_particle (qmax_org_fine_particle was calculated by Lin (lyu@bgc-jena.mpg.de) per NWRB category available from TAXNWRB.soilgrid.3600.1800.nc)
#  only update lon_sg/lat_sg but use the data objects alread read above (all_cruncep section)
#
inNWRBcatFile <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/TAXNWRB.soilgrid.3600.1800.nc"

ncid         <- nc_open(inNWRBcatFile)  # grid with 3600x1800 resolution
lon_sg       <- ncvar_get(ncid,'longitude')
lat_sg       <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create empty vectors for site-data
qmaxorg_site <- numeric(sppSiteNum)

## loop over the two sites: get NWRB cat data from nearest location & get according Qmax_org
for (s in 1:sppSiteNum) {
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
## prepare the output file: spp_gfdb_sites* 1901-2015 and 1984-2013
##  

# output file
siteParam <- 27
siteinfo  <- array(NA,c(sppSiteNum,siteParam))

# loop over gfdb spp selected sites
i <- 1
for (spp_site in spp_site_list) {

    siteinfo[i,1]   <-  paste0(spp_site,"_01")
    siteinfo[i,2]   <-  lon_orig[i]
    siteinfo[i,3]   <-  lat_orig[i]
    siteinfo[i,4]   <-  pft_spp[i]
    siteinfo[i,5]   <-  1901
    siteinfo[i,6]   <-  2015
    siteinfo[i,7]   <-  lon_orig[i]
    siteinfo[i,8]   <-  clay[i]
    siteinfo[i,9]   <-  silt[i]
    siteinfo[i,10]  <-  sand[i]
    siteinfo[i,11]  <-  awc[i]
    siteinfo[i,12]  <-  bulk_density[i]
    siteinfo[i,13]  <-  ph[i]
    siteinfo[i,14]  <-  taxousda[i]
    siteinfo[i,15]  <-  taxnwrb[i]
    siteinfo[i,16]  <-  glim_class_site[i]
    siteinfo[i,17]  <-  (-9999) #LAI
    siteinfo[i,18]  <-  (-9999) #Nleaf
    siteinfo[i,19]  <-  (-9999) #SLA
    siteinfo[i,20]  <-  (-9999) #Height
    siteinfo[i,21]  <-  1890
    siteinfo[i,22]  <-  soilP_depth_site[i]
    siteinfo[i,23]  <-  soilP_labile_site[i]
    siteinfo[i,24]  <-  soilP_slow_site[i]
    siteinfo[i,25]  <-  soilP_occlud_site[i]
    siteinfo[i,26]  <-  soilP_primary_site[i]
    siteinfo[i,27]  <-  qmaxorg_site[i]

  i<-i+1
}


## output
# time stamp
date <- Sys.Date()
# header
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')

## CRUNCEP SPP
# 1901-2015
write.table(siteinfo,file=paste0("crun_spp_sites_1901-2015_list",date,".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
# 1984-2013
write.table(siteinfo,file=paste0("crun_spp_sites_1984-2013_list",date,".dat"),quote=F,row.names=F,col.names=T)

## SPP_1685
i <- 1
for (spp_site in spp_site_list) {
  siteinfo[i,1]   <-  spp_site  # site_ID without an appended "_01"
  i<-i+1
}
# 1984-2013
write.table(siteinfo,file=paste0("spp_1685_sites_1984-2013_list",date,".dat"),quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1901
siteinfo[,6]<-2015
# 1901-2015
write.table(siteinfo,file=paste0("spp_1685_sites_1901-2015_list",date,".dat"),quote=F,row.names=F,col.names=T)



print(" ")
print("Script has finished.")
print(paste0("The new crun_all/cue/spp_*",date,".dat & spp_1685_sites_*",date,".dat files with ",siteParam," Parameters have been created."))


## end ##
