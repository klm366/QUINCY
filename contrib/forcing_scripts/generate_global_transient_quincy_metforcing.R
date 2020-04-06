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
## generate > transient < met forcing for a set of sites GLOBAL distributed across the globe (similar to CRUNCEP site-set)
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_global_transient_quincy_metforcing.R
## 
## note: allsiteinfo is ./generate_global_all_site_list.dat  (columns: Site-ID lon lat PFT start end)
## 
## output may be generated per default in:
##  ./new_gobal_transient_met_forcing/
##  output file contains "_t_" as identifier for transient forcing
## 
#########################################################################################


## library
library(ncdf4)

## create out folder & set in/out path
localOutDir <- "./new_gobal_transient_met_forcing"
dir.create(localOutDir, showWarnings=T, recursive=T)
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/GLOBAL_siteset_data/"
ipath2 <- "/Weather/"              # relative path within each site folder
infile_prefix <- 'weather_dump_'
infile_ext    <- '.nc'
opath  <- paste0(localOutDir, "/")


## forcing period
startyear <- 1901
endyear   <- 2018

## site info
allsiteinfo <- read.table("generate_global_all_site_list.dat", stringsAsFactor=F, header=T, sep="")  # (columns: Site-ID lon lat PFT start end)
allsitename <- allsiteinfo[,1]
allsitelon  <- allsiteinfo[,2]
allsitelat  <- allsiteinfo[,3]
site_list   <- allsitename


doy=ceiling(seq(1/48,365,1/48))
hod=rep(seq(0,23.5,0.5),365)
dpm<-c(31,28,31,30,31,30,31,31,30,31,30,31)

header1<-c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
           'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2<-c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

co2y<-read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

## loop over site list
for (s in 1:length(site_list)) {
  print(paste0("site: ", site_list[s], "  number ", s, " of ", length(site_list)))
  # check if site data exists by checking if the year 1994 forcing file exists
  stat3   <- file.exists(paste0(ipath1, site_list[s], ipath2, infile_prefix, 1994, infile_ext))
  # name of output file contains "_t_" as identifier for transient forcing
  outname <- paste0(opath, site_list[s], "_t_", startyear, "-", endyear, ".dat")
  # if site data exists
  if(stat3) {
    if (file.exists(outname)) file.remove(outname) 
    # write header
    write.table(t(cbind(header1,header2)),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ') 
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      # meteo data
      fid<-nc_open(paste0(ipath1,site_list[s],ipath2,infile_prefix,year,infile_ext))
       swdown<-ncvar_get(fid,"SWdown")
       lwdown<-ncvar_get(fid,"LWdown")
       tair<-ncvar_get(fid,"Tair")
       rain<-ncvar_get(fid,"Rainf")*86400
       snow<-ncvar_get(fid,"Snowf")*86400
       qair<-ncvar_get(fid,"Qair")*1000
       press<-ncvar_get(fid,"PSurf")/100
       wind<-ncvar_get(fid,"Wind_N")
      nc_close(fid)
      # N deposition data
      # ndep_YYYY.nc        available 1850-1995
      # ndep_RCP8.5_YYYY.nc available 1996-2100  # also other RCP available
      if (year < 1996){
        nDepFilePref <- "ndep_"
      }else{
        nDepFilePref <- "ndep_RCP8.5_"
      }
      fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/", nDepFilePref, year, ".nc"))
       lonndep<-ncvar_get(fid,'lon')
       latndep<-ncvar_get(fid,'lat')
       lonsite<-allsitelon[site_list[s]==allsitename]
       latsite<-allsitelat[site_list[s]==allsitename]
       latx<-which.min(abs(latndep-latsite))
       lonx<-which.min(abs(lonndep-lonsite))
       nhx_dummy<-ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
       noy_dummy<-ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
       noy=nhx<-numeric(length(swdown))
       nhx[1:(dpm[1]*48)]<-rep(nhx_dummy[1],dpm[1]*48)
       noy[1:(dpm[1]*48)]<-rep(noy_dummy[1],dpm[1]*48)
       for (m in 2:12) {
         nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(nhx_dummy[m],dpm[m]*48)
         noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(noy_dummy[m],dpm[m]*48)
       }
      nc_close(fid)
      # P deposition data
      fid<-nc_open(paste0("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc"))
       lonpdep<-ncvar_get(fid,'lon')
       latpdep<-ncvar_get(fid,'lat')
       lonsite<-allsitelon[site_list[s]==allsitename]
       latsite<-allsitelat[site_list[s]==allsitename]
       latx<-which.min(abs(latpdep-latsite))
       lonx<-which.min(abs(lonpdep-lonsite))
       p_dummy<-ncvar_get(fid,'pdep')[lonx,latx]/365
       pdep<-numeric(length(swdown))
       pdep[]<-p_dummy
      nc_close(fid)
      #CO2
      co2<-rep(co2y[co2y[,1]==year,2],length(swdown))
      co2_dC13<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],length(swdown))
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      co2_DC14<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],length(swdown))
      year_out<-rep(year,length(swdown))
      # combine into one object
      forcing<-cbind(year_out,doy,hod,swdown,lwdown,tair,qair,press,rain,snow,wind,co2,co2_dC13,co2_DC14,nhx,noy,pdep)
      # write to file, append=T
      write.table(format(round(forcing,digits=4),nsmall=1),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ')
    }
  } else {
    stat2<-unlist(strsplit(stat,split=' '))[11]
    filestat<-unlist(strsplit(stat2,split="/"))
    linktarget<-filestat[length(filestat)-1]
    refname=paste0(linktarget,"_",startyear,"-",endyear,'.dat')
       print(c('linking',outname,refname))
       system(paste0("ln -s",refname,outname))   
  }
}


## fin ##
