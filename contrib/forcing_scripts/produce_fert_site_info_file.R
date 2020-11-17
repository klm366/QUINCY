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
## generate site_list.dat files for the FERT site set for QUINCY
##

## info:
##  * the FERT site set is used for two different experiments: a) N fertilization, and b) 15N pulses (site COR or CR2)
##  * a) the site info file is created for "Site-ID"-CON sites, with the "CON" tag for control
##    * these sites are the "base set" of sites 
##    * manually added are/will be (by Silvia) copies of each of these sites, with tags like "-HIG" "-LOW" "-ADD" etc.
##    * again, this script does only screate a site info file with ""-CON" basic sites
##  * b) the COR (Cornell Uni) is used for the experiment
##    * the 15N pulses are not added here (nor are they added with the generate_fert_transient_quincy_forcing.R)
##    * the idea is to add it at runtime in the QUINCY model, probably hardcoded specific to a site
##    * hence the same COR site is present in this site-set two times with a different name each ("-CON" & "-15N" tags)
##    * in addition, new forcing data for the COR site was available with April 2020, hence, 
##      another two COR sites were added to this site set named CR2 for "COR version 2", with forcing to yr 2018
##    * for the CR2 site two versions are created, same as for COR ("-CON" & "-15N" tags)

## library
library(ncdf4)

#########################################################################################
## run this script by:
## cd contrib/forcing_scripts/
## Rscript --vanilla ./produce_fert_site_info_file.R


#########################################################################################
## output 
## the output files will be written in this folder


#########################################################################################
## basic site list
metsites_fert     <- read.table("generate_fert_siteinfo_noheader.dat", stringsAsFactors=F, header=F, sep="")
fert_site_id_list <- metsites_fert[,1]     # Site-ID (the "COR" & "CR2" sites are present two times; see above info)
allSiteNameD      <- metsites_fert[,23]    # detailed sites names (last column)

# define start and last year
# last year is 2015 by default, but 2018 for the CR2 site(s)
startYear <- rep(1901, length(fert_site_id_list))
lastYear  <- rep(2015, length(fert_site_id_list))
lastYear[which(fert_site_id_list == "CR2")] <- 2018


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


#########################################################################################
## preparing & creating the output file
## 
## 
## new data is read from files / generated here
##
##

## get number of sites (see file input above; no column header)
fertSiteNum <- length(fert_site_id_list)

## get the lon/lat of all sites (but VES) from the metsites_fert object (with list of siteID in metsites_fert), file siteinfo_noheader.dat
#  new values for lon_orig & lat_orig  (above they contain lon/lat of all sites in forestsite_list from info_forest.csv)
lon_orig <- numeric(fertSiteNum)
lat_orig <- numeric(fertSiteNum)
for (s in 1:length(fert_site_id_list)) {
  lon_orig[s] <- metsites_fert[s,2]
  lat_orig[s] <- metsites_fert[s,3]
}


## update the lon/lat values of the soilgrid used below for deriving soil data
#  only update lon_sg/lat_sg but use the data objects alread read above
#
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/"

ncid   <- nc_open(paste0(ipath,"SNDPPT.soilgrid.3600.1800.7.nc"))
lon_sg <- ncvar_get(ncid,'longitude')
lat_sg <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create new (empty) vectors/arrays for output
depth_begin  <- c(0.0,0.05,0.15,0.3,0.6,1.0,2.0)
depth_end    <- c(0.05,0.15,0.3,0.6,1.0,2.0,100)
depth_sl     <- depth_end - depth_begin
awc          <- numeric(fertSiteNum)
clay         <- numeric(fertSiteNum)
silt         <- numeric(fertSiteNum)
sand         <- numeric(fertSiteNum)
taxnwrb      <- numeric(fertSiteNum)
taxousda     <- numeric(fertSiteNum)
depth        <- numeric(fertSiteNum)
ph           <- numeric(fertSiteNum)
bulk_density <- numeric(fertSiteNum)
awc_sl       <- array(NA,c(fertSiteNum,7))

## loop over the sites and get soil data from above read netcdf files
for (s in 1:fertSiteNum) {
  lonidx <- which.min(abs(lon_sg-lon_orig[s]))
  latidx <- which.min(abs(lat_sg-lat_orig[s]))
  
  bulk_density[s] <- mean(bulk_density_sg[lonidx,latidx,])
  depth[s]        <- depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]     <- taxousda_sg[lonidx,latidx]
  taxnwrb[s]      <- taxnwrb_sg[lonidx,latidx]
  awc_sl[s,]      <- (fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  awc[s]          <- weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
  
  ph[s]           <- mean(ph_sg[lonidx,latidx,])/10
  sand[s]         <- weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  silt[s]         <- weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  clay[s]         <- 1.0 - sand[s] - silt[s]
}

## PFT
## set PFTs according to the OCN model output generated by Sönke (object metsites_fert from generate_fert_siteinfo_noheader.dat)
## this site set does contain only PFT 4 & 5 sites
## info from Sönke: PFT is 4 if column 15 has a 1; PFT is 5 if column 13 has a 1
## QUINCY model:    PFT 4: TeBS  /                 PFT 5: TeNE
# sites names are saved in: fert_site_id_list (sites currently are: HAR SDA NRR NRL ASE FAR LLA LOV ASA SKB MAS MIB MIC MID HRL OSS COR COR)
# create a PFT vector with PFT 4 TeBS as default
pft_fert <- c(rep("TeBS", length(fert_site_id_list)))
# loop over sites and change PFT to 5 where necessary
for (s in 1:length(fert_site_id_list)){
  # PFT 4 TeBS if column 15 = 1
    # default
  # PFT 5 TeNE if column 13 = 1
  if (metsites_fert[s,13] == 1) pft_fert[s] = "TeNE"
}

## read lithology map
#  this is the same input file as above, but repeating it here again for improved modularity
#
ipath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/"

ncid            <- nc_open(paste0(ipath,"GLim.720.360.nc"))
lon_glim        <- ncvar_get(ncid,'lon')
lat_glim        <- ncvar_get(ncid,'lat')
glim_class      <- ncvar_get(ncid,'GLim')
glim_class_site <- numeric(fertSiteNum)
nc_close(ncid)

## loop over sites: get soil data
for (s in 1:fertSiteNum) {
  lonidx             <- which.min(abs(lon_glim-lon_orig[s]))
  latidx             <- which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s] <- glim_class[lonidx,latidx]
}


## read site-specific data on Phosphorous
#  only update lon_sg/lat_sg but use the data objects alread read above
#
inPhosPath <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/Phosphorous/v2014_06/Data/"

ncid             <- nc_open(paste0(inPhosPath,"Labile_Inorganic_P.nc"))
lon_sg           <- ncvar_get(ncid,'longitude')
lat_sg           <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create empty vectors for site-data
soilP_depth_site   <- numeric(fertSiteNum)
soilP_labile_site  <- numeric(fertSiteNum)
soilP_slow_site    <- numeric(fertSiteNum)
soilP_occlud_site  <- numeric(fertSiteNum)
soilP_primary_site <- numeric(fertSiteNum)

## loop over the two sites: get soil P data from nearest location
for (s in 1:fertSiteNum) {
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
#  only update lon_sg/lat_sg but use the data objects already read above
#
inNWRBcatFile <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/TAXNWRB.soilgrid.3600.1800.nc"

ncid         <- nc_open(inNWRBcatFile)  # grid with 3600x1800 resolution
lon_sg       <- ncvar_get(ncid,'longitude')
lat_sg       <- ncvar_get(ncid,'latitude')
nc_close(ncid)

## create empty vectors for site-data
qmaxorg_site <- numeric(fertSiteNum)

## loop over the sites: get NWRB cat data from nearest location & get according Qmax_org
for (s in 1:fertSiteNum) {
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
## prepare the output file: fert_sites* 1901-2015
##  

# output file
siteParam <- 27
siteinfo  <- array(NA,c(fertSiteNum,siteParam))

# loop over all FERT sites
i <- 1
for (fert_site in fert_site_id_list) {

    siteinfo[i,1]   <-  paste0(fert_site,"-CON")
    # site name special case: COR site manipulated for 15N deposition pulse experiment
    if (fert_site == "COR" & allSiteNameD[i] == "15N_experiment")      siteinfo[i,1] <- "COR-15N"
    if (fert_site == "CR2" & allSiteNameD[i] == "15N_experiment_2018") siteinfo[i,1] <- "CR2-15N"
    siteinfo[i,2]   <-  lon_orig[i]
    siteinfo[i,3]   <-  lat_orig[i]
    siteinfo[i,4]   <-  pft_fert[i]
    siteinfo[i,5]   <-  1901
    siteinfo[i,6]   <-  2015
    siteinfo[i,7]   <-  lon_orig[i]
    siteinfo[i,8]   <-  round(clay[i], 5)
    siteinfo[i,9]   <-  round(silt[i], 5)
    siteinfo[i,10]  <-  round(sand[i], 5)
    siteinfo[i,11]  <-  round(awc[i], 5)
    siteinfo[i,12]  <-  round(bulk_density[i], 5)
    siteinfo[i,13]  <-  round(ph[i], 5)
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
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end',
                      'lon_gmt','clay','silt','sand','awc','bd','ph',
                      'taxusda','taxnwrb','lith_glim',
                      'LAI','Nleaf','SLA','Height','PlantYear',
                      'soilP_depth','soilP_labile','soilP_slow','soilP_occluded','soilP_primary',
                      'Qmax_org_fp')
# 1901-2015 for all sites but CR2
# 1901-2018 for CR2
write.table(siteinfo,file=paste0("fert_sites_1901-2018_list",Sys.Date(),".dat"), quote=F, row.names=F, col.names=T)


print(" ")
print(paste0(" Script has finished. The new fert_sites_1901-2018_listDATE.dat files with ",siteParam," Parameters have been created."))
print(" ")
print(" Please note: this basic version of a site info file may be manually extended to more sites, by just multiplying some sites and give it diff names.")
print("   For more details, see the file already in place: input/point/FERT/fert_sites_1901-2018_list.dat")

## end ##
