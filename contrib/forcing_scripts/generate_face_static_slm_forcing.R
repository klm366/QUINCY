#!/usr/bin/env Rscript

##    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
##    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
##
##    This file is part of the QUINCY model.
##    The QUINCY model is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
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
## generate > static < met forcing for FACE sites
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_face_static_slm_forcing.R
## 
## note: allsiteinfo is ./generate_face_site_list.dat 
## 
## output may be generated per default in:
##  ./new_face_static_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)




## create out folder & set in/out path
dir.create("./new_face_static_met_forcing", showWarnings=T, recursive=T)
ipath1='/Net/Groups/BGI/model/OCN/input/point/FACE/'
infile_ext='_FACE.nc'
opath="./new_face_static_met_forcing/"
#opath="/Net/Groups/BGI/data/Model/SLM1/input/point/FACE/met_forcing/"

## read site info file
allsiteinfo<-read.table("generate_face_site_list.dat", stringsAsFactor=F, header=T)  # file is originally from bgi/people/szaehle/Projects/QUINCY/GenerateFACEForcing_2017/site_list.dat
allsitename<-allsiteinfo[,1]
allsitestartyear<-allsiteinfo[,2]
allsiteendyear<-allsiteinfo[,3]
allsitelevel<-allsiteinfo[,5]
allsite_newname<-allsiteinfo[,6]
allsite_ed13c<-allsiteinfo[,7]
allsite_ed14c<-allsiteinfo[,8]
allsite_ndep<-allsiteinfo[,9]
site_list<-allsitename


doy=ceiling(seq(1/48,365,1/48))
hod=rep(seq(0,23.5,0.5),365)
dpm<-c(31,28,31,30,31,30,31,31,30,31,30,31)

header1<-c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
           'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2<-c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

co2y<-read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

slm_site_list<-array(0,c(length(site_list),17))

for (s in 1:length(site_list)) {
  if(allsitelevel[s]=='AMB') level=3
  if(allsitelevel[s]=='ELE') level=4
  co2daily<-read.table(paste(ipath1,'co2_daily_',site_list[s],'.dat',sep=''))
#for (s in 1:1) {
print(site_list[s])
  startyear<-allsitestartyear[allsite_newname==allsite_newname[s]]
  endyear<-allsiteendyear[allsite_newname==allsite_newname[s]]
  nyear<-endyear-startyear+1
  stat3<-file.exists(paste(ipath1,site_list[s],infile_ext,sep=''))
  outname=paste(opath,allsite_newname[s],"_",startyear,"-",endyear,'.dat',sep='')
  if(stat3) {
    if (file.exists(outname)) file.remove(outname) 
    # write header
    write.table(t(cbind(header1,header2)),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ') 
      # meteo data
      fid<-nc_open(paste(ipath1,site_list[s],infile_ext,sep=''))
       lonsite<-ncvar_get(fid,"x")
       latsite<-ncvar_get(fid,"y")
       swdown<-ncvar_get(fid,"SWdown")
       lwdown<-ncvar_get(fid,"LWdown")
       tair<-ncvar_get(fid,"Tair")
       rain<-ncvar_get(fid,"Rainf")*86400
       snow<-ncvar_get(fid,"Snowf")*86400
       qair<-ncvar_get(fid,"Qair")*1000
       press<-ncvar_get(fid,"PSurf")/100
       wind<-ncvar_get(fid,"Wind_N")
       # something is odd with the wind forcing of Duke for the first two years
       tsel<-seq(3,10,1)*17520
       if(site_list[s]=='DUKE') {
          for (i in 1:35088) wind[i]<-mean(wind[i-1+tsel]) 
       }
       wind[wind<0.1]<-0.1
      nc_close(fid)
      # P deposition data
      fid<-nc_open(paste("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc",sep=''))
       lonpdep<-as.numeric(ncvar_get(fid,'lon'))
       latpdep<-as.numeric(ncvar_get(fid,'lat'))
       latx<-which.min(abs(latpdep-latsite))
       lonx<-which.min(abs(lonpdep-lonsite))
       p_dummy<-ncvar_get(fid,'pdep')[lonx,latx]/365
       pdep<-numeric(length(swdown))
       pdep[]<-p_dummy
      nc_close(fid)
      # N deposition data
    year_out=co2=co2_dC13=co2_DC14=noy=nhx<-numeric(length(swdown))
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      fid<-nc_open(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc",sep=''))
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
      #CO2
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      if(year>=startyear) {
         for (d in 1:365) { 
            co2[1:48+(d-1)*48+365*48*(year-startyear)]<-co2daily[co2daily[,1]==year&co2daily[,2]==d,level]
            if (level==4&(co2daily[co2daily[,1]==year&co2daily[,2]==d,3]!=co2daily[co2daily[,1]==year&co2daily[,2]==d,4])) {
               co2_dC13[1:48+(d-1)*48+365*48*(year-startyear)]<-(co2daily[co2daily[,1]==year&co2daily[,2]==d,3]*co2y_dC13[(co2y_dC13[,1]-0.5)==year,2]+
                                         (co2daily[co2daily[,1]==year&co2daily[,2]==d,4]-co2daily[co2daily[,1]==year&co2daily[,2]==d,3])*allsite_ed13c[s])/co2daily[co2daily[,1]==year&co2daily[,2]==d,4]
               co2_DC14[1:48+(d-1)*48+365*48*(year-startyear)]<-(co2daily[co2daily[,1]==year&co2daily[,2]==d,3]*co2y_DC14[(co2y_dC13[,1]-0.5)==year,2]+
                                         (co2daily[co2daily[,1]==year&co2daily[,2]==d,4]-co2daily[co2daily[,1]==year&co2daily[,2]==d,3])*allsite_ed14c[s])/co2daily[co2daily[,1]==year&co2daily[,2]==d,4]
            } else {
               co2_dC13[1:48+(d-1)*48+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],48)
               co2_DC14[1:48+(d-1)*48+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],48)
            }
         }
      }
      else {
         co2[1:17520+365*48*(year-startyear)]<-rep(co2y[co2y[,1]==year,2],17520)
         co2_dC13[1:17520+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],17520)
         if(latsite>30) C14lat=2
         if(latsite>-30&latsite<=30) C14lat=3
         if(latsite<=-30) C14lat=4
         co2_DC14[1:17520+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],17520)
      }
      year_out[1:17520+365*48*(year-startyear)]<-rep(year,17520)
      }
      doy_out<-(rep(doy,nyear))
      hod_out<-(rep(hod,nyear))
      nstep<-length(doy_out)
      forcing<-cbind(year_out[1:nstep],doy_out,hod_out,swdown[1:nstep],lwdown[1:nstep],tair[1:nstep],qair[1:nstep],
                     press[1:nstep],rain[1:nstep],snow[1:nstep],wind[1:nstep],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                     nhx[1:nstep],noy[1:nstep],pdep[1:nstep])
      write.table(format(round(forcing,digits=4),nsmall=1),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ')
  } else {
       print('site not available')
  }
 
}
