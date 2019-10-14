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
## generate > transient < met forcing for Cruncep v7 SPP (& SPP_1685) sites using weatherdump v2
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_cruncep_spp_transient_slm_forcing.R
## 
## note: allsiteinfo is ./generate_cruncep_spp_siteinfo_noheader.dat 
## 
## output may be generated per default in:
##  ./new_cruncep_spp_transient_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)

## create out folder & set in/out path
dir.create("./new_cruncep_spp_transient_met_forcing", showWarnings=T, recursive=T)
ipath1="/Net/Groups/BGI/people/szaehle/Projects/QUINCY/GenerateCRUForcing_2017/Results/"
ipath2="/Weather/"              # relative path within each site folder
infile_prefix='weather_dump_'
infile_ext='.nc'
opath="./new_cruncep_spp_transient_met_forcing/"
#opath='/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv7/met_forcing/'

## forcing period
startyear <- 1901
endyear   <- 2015

## site info
allsiteinfo <- read.table("generate_cruncep_spp_siteinfo_noheader.dat", stringsAsFactor=F)
allsitelon  <- allsiteinfo[,2]
allsitelat  <- allsiteinfo[,3]
allsitename <- allsiteinfo[,1]
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

## N deposition per site [kg N/ha]
# using this number allows for scaling the N deposition, i.e., consider site history / pre-industrial values
# these numbers come from literature:
#  BBR and LUE: Bünemann, E. K., Augstburger, S. & Frossard, E. (2016) Dominance of either physicochemical or biological phosphorus cycling processes 
#               in temperate forest soils of contrasting phosphate availability. Soil Biology and Biochemistry 101, 85-95, doi:10.1016/j.soilbio.2016.07.005
#  VES: Forstliches Umweltmonitoring in Thüringen (Gotha, 2001)
#  CON: calculated from raw data; data source: Jaane Krüger, personal communication with Lin Yu
#  MTF: no literature value available; calculated as the mean value of CON, BBR and LUE
# BBR & LUE & CON are about 2014 values and VES is the value from the year 2000
# BBR LUE VES MTF CON
ndep_site<-c(8, 13, 20.65, 9.7, 8.19)

for (s in 1:length(site_list)) {
print(paste0("site: ", site_list[s], "  number ", s, " of ", length(site_list)))
  orig_climate_file_name<-paste("/Net/Groups/BGI/model/OCN/input/point/CRUNCEP/",site_list[s],"/climate_1994.nc",sep='')
  #if(file.exists(orig_climate_file_name)) {
  #   stat<-system(paste("ls -l ",orig_climate_file_name,sep=''),intern=T) 
  #   sizestat<-as.numeric(unlist(strsplit(stat,split=' '))[5])
  #} else {sitestat=0}
  stat3<-file.exists(paste(ipath1,site_list[s],ipath2,infile_prefix,1994,infile_ext,sep=''))
  outname=paste(opath,site_list[s],"_",startyear,"-",endyear,'.dat',sep='')
  if(stat3) {

    # calculate N deposition scalar to correct N deposition to local observed value
    fid<-nc_open(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_1861.nc",sep=''))
    fid2<-nc_open(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_2005.nc",sep=''))
     lonndep<-ncvar_get(fid,'lon')
     latndep<-ncvar_get(fid,'lat')
     lonsite<-allsitelon[site_list[s]==allsitename]
     latsite<-allsitelat[site_list[s]==allsitename]
     latx<-which.min(abs(latndep-latsite))
     lonx<-which.min(abs(lonndep-lonsite))
     n_preindus<-mean(ncvar_get(fid,'NHx_deposition')[lonx,latx,]+ncvar_get(fid,'NOy_deposition')[lonx,latx,])*86400*365*10000
     nhx_preindus<-ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
     noy_preindus<-ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
     n_present<-mean(ncvar_get(fid2,'NHx_deposition')[lonx,latx,]+ncvar_get(fid2,'NOy_deposition')[lonx,latx,])*86400*365*10000
    nc_close(fid2)
    nc_close(fid)

    # overwrite existing file
    if (file.exists(outname)) file.remove(outname) 
    # write header
    write.table(t(cbind(header1,header2)),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ') 
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      # meteo data
      fid<-nc_open(paste(ipath1,site_list[s],ipath2,infile_prefix,year,infile_ext,sep=''))
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
      fid<-nc_open(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc",sep=''))
       lonndep<-ncvar_get(fid,'lon')
       latndep<-ncvar_get(fid,'lat')
       lonsite<-allsitelon[site_list[s]==allsitename]
       latsite<-allsitelat[site_list[s]==allsitename]
       latx<-which.min(abs(latndep-latsite))
       lonx<-which.min(abs(lonndep-lonsite))
       scal<-(ndep_site[s]-n_preindus)/(n_present-n_preindus)
       nhx_dummy<-scal*(ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000-nhx_preindus)+nhx_preindus
       noy_dummy<-scal*(ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000-noy_preindus)+noy_preindus
       n_year<-mean(nhx_dummy+noy_dummy)*365/100
       print(c(scal,n_year,ndep_site[s],n_present))
       noy=nhx<-numeric(length(swdown))
       nhx[1:(dpm[1]*48)]<-rep(nhx_dummy[1],dpm[1]*48)
       noy[1:(dpm[1]*48)]<-rep(noy_dummy[1],dpm[1]*48)
       for (m in 2:12) {
         nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(nhx_dummy[m],dpm[m]*48)
         noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(noy_dummy[m],dpm[m]*48)
       }
      nc_close(fid)
      # P deposition data
      fid<-nc_open(paste("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc",sep=''))
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
      forcing<-cbind(year_out,doy,hod,swdown,lwdown,tair,qair,press,rain,snow,wind,co2,co2_dC13,co2_DC14,nhx,noy,pdep)
      write.table(format(round(forcing,digits=4),nsmall=1),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ')
    }
  } else {
    stat2<-unlist(strsplit(stat,split=' '))[11]
    filestat<-unlist(strsplit(stat2,split="/"))
    linktarget<-filestat[length(filestat)-1]
    refname=paste(linktarget,"_",startyear,"-",endyear,'.dat',sep='')
       print(c('linking',outname,refname))
       system(paste("ln -s",refname,outname))   
  }
 
}
