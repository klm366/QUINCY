#!/usr/bin/env Rscript

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

#########################################################################################
##
## generate > static < met forcing for four finish sites
##
## the original forcing input does contain 5 leap years (1996 2000 2004 2008 2012)
##  in the below script all dates "February 29th" are neglected, i.e., if time = YYYY0229.x
## 
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_confin_static_quincy_forcing.R
## 
## note: * allsiteinfo is ./generate_confin_site_list.dat  
##         (manually created from BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/site_descriptions.txt)

# 2do !! <<<<<<<<<<<<<<<<<<<<< check & do the site info file creation        <<<<<<<<<<<<<  !!!!

##       * also the actual site info file (confin_sites_1996-2014_list.dat) was created manually
##         (from files: jsbach_def_4sites.nc ? and site_descriptions.txt in BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/)
##         using the fluxnet FI-Hyy (from FLUXNET2 site-set) as basis that is modified
##       * this site info file needs (manual) modification when a new dataset with more recent data is processed
## 
## output may be generated per default in:
##  ./new_confin_static_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)


## create out folder & set in/out path
dir.create("./new_confin_static_met_forcing", showWarnings=T, recursive=T)
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/"
opath  <- "./new_confin_static_met_forcing/"

## define & load climate-forcing input file
#  contains the data for all four sites
infile <- "climate_4sites_1996-2014_with_leap_years.nc"
load(paste0(ipath1, infile))

## define first and last year of dataset
firstYear <- 1996
lastYear  <- 2014
leapYears <- c(1996, 2000, 2004, 2008, 2012)

## message
print(paste0("Generating met forcing for four finish conifers sites (CONFIN) for : ", firstYear, "to ",lastYear))

## read site info file (manually created from BGI/people/szaehle/Projects/QUINCY/GenerateMANP_Forcing_2018/site_info.dat)
allsiteinfo       <- read.table("generate_confin_site_list.dat", stringsAsFactor=F, header=T)
allsitestartyear  <- allsiteinfo[,5]
allsiteendyear    <- allsiteinfo[,6]
allsitename       <- allsiteinfo[,1]
site_list         <- allsitename

doy=ceiling(seq(1/48,365,1/48))
hod=rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)  # days per month

forcVar <- 17  # number of forcing variables including time (see header1 & 2)
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
           'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

co2y      <- read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

slm_site_list <- array(0, c(length(site_list), forcVar))

for (s in 1:length(site_list)) {
#for (s in 1:1) {
print(site_list[s])
  startyear<-allsitestartyear[allsitename==site_list[s]]
  endyear<-allsiteendyear[allsitename==site_list[s]]
  nyear<-endyear-startyear+1
  outname=paste(opath,site_list[s],"_",startyear,"-",endyear,'.dat',sep='')
  #if(stat3) {
    if (file.exists(outname)) file.remove(outname) 
    # write header
    write.table(t(cbind(header1,header2)),file=outname,append=TRUE,
                quote=FALSE,col.names=FALSE,row.names=FALSE,sep=' ') 
      # meteo data
       lonsite<-allsiteinfo[allsitename==site_list[s],2]
       latsite<-allsiteinfo[allsitename==site_list[s],3]
       swdown<-df[,4]
       lwdown<-df[,5]
       tair<-df[,6]
       rain<-df[,9]
       snow<-df[,10]
       qair<-df[,7]
       press<-df[,8]
       wind<-df[,11]
       wind[wind<0.1]<-0.1
      # P deposition data
      fid<-nc_open("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc")
       lonpdep<-as.numeric(ncvar_get(fid,'lon'))
       latpdep<-as.numeric(ncvar_get(fid,'lat'))
       latx<-which.min(abs(latpdep-latsite))
       lonx<-which.min(abs(lonpdep-lonsite))
       p_dummy<-ncvar_get(fid,'pdep')[lonx,latx]/365
       pdep<-numeric(length(swdown))
       pdep[]<-p_dummy
      nc_close(fid)

    ## add data looping over the years
    year_out=co2=co2_dC13=co2_DC14=noy=nhx<-numeric(length(swdown))
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      ## N deposition
      #  using RCP 8.5 data throughout
       fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_RCP8.5_",year,".nc"))
       lonndep<-ncvar_get(fid,'lon')
       latndep<-ncvar_get(fid,'lat')
       latx<-which.min(abs(latndep-as.numeric(latsite)))
       lonx<-which.min(abs(lonndep-as.numeric(lonsite)))
       nhx_dummy<-ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
       noy_dummy<-ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
       nhx[1:(dpm[1]*48)+365*48*(year-startyear)]<-rep(nhx_dummy[1],dpm[1]*48)
       noy[1:(dpm[1]*48)+365*48*(year-startyear)]<-rep(noy_dummy[1],dpm[1]*48)
       for (m in 2:12) {
         nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(nhx_dummy[m],dpm[m]*48)
         noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(noy_dummy[m],dpm[m]*48)
       }
      nc_close(fid)
      ## CO2 
      #  use global mean values from GCP_co2_global.dat for all values before 2014,
      #  because the Majadas field sensor was not well calibrated before 2014 and, hence, does not give reasonable values
      if(year<2014) {
        co2[1:17520+365*48*(year-startyear)]<-rep(co2y[co2y[,1]==year,2],17520)
      } else {
        co2[1:17520+365*48*(year-startyear)]<-rep(co2y[co2y[,1]==2016,2],17520)+2*(year-2016)
      }
      if(year<2014) {
        co2_dC13[1:17520+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],17520)
      } else {
        co2_dC13[1:17520+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==2015,2],17520)-0.02*(year-2015)
      }
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      if(year<2014) {
        co2_DC14[1:17520+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],17520)
      } else {
        co2_DC14[1:17520+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==2015,C14lat],17520)-5*(year-2015)
      }
      year_out[1:17520+365*48*(year-startyear)]<-rep(year,17520)
      }
      doy_out<-(rep(doy,nyear))
      hod_out<-(rep(hod,nyear))
      nstep<-length(doy_out)
      
      # combine individual variables into one array for output
      forcing<-cbind(year_out[1:nstep],doy_out,hod_out,swdown[1:nstep],lwdown[1:nstep],tair[1:nstep],qair[1:nstep],
                     press[1:nstep],rain[1:nstep],snow[1:nstep],wind[1:nstep],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                     nhx[1:nstep],noy[1:nstep],pdep[1:nstep])
                     
      # write output append yearly data (header created above)
      write.table(format(round(forcing,digits=4),nsmall=1), file=outname, append=TRUE, quote=FALSE,
                  col.names=FALSE, row.names=FALSE, sep=" ")

}


## message
print("Done !")


## end ##
