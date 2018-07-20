#!/usr/bin/env Rscript
library(ncdf)
load("ES_LMa_2009-2017_gf.Rdata")

opath="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/met_forcing/"

allsiteinfo<-read.table('site_info.dat',stringsAsFactor=F,header=T)
allsitestartyear<-allsiteinfo[,5]
allsiteendyear<-allsiteinfo[,6]
allsitename<-allsiteinfo[,1]
site_list<-allsitename

doy=ceiling(seq(1/48,365,1/48))
hod=rep(seq(0,23.5,0.5),365)
dpm<-c(31,28,31,30,31,30,31,31,30,31,30,31)

header1<-c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
           'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2<-c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

header_site_list<-c('Site-ID','lon','lat','PFT','start','end','lon_gmt','clay','silt','sand','awc','pH','LAI','Nleaf','SLA','Height','PlantYear')

co2y<-read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14<-read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

slm_site_list<-array(0,c(length(site_list),17))

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
    write.table(t(cbind(header1,header2)),file=outname,append=TRUE,quote=FALSE,col.names=FALSE,
              row.names=FALSE,sep=' ') 
      # meteo data
  #    fid<-open.ncdf(paste(ipath1,site_list[s],'.',startyear,'-',endyear,'.',infile_ext,sep=''))
       lonsite<-allsiteinfo[allsitename==site_list[s],2]
       latsite<-allsiteinfo[allsitename==site_list[s],3]
       swdown<-df[,4]
       lwdown<-df[,5]
       tair<-df[,6]
       rain<-df[,9]*48*48
       snow<-df[,10]
       qair<-df[,7]
       press<-df[,8]
       wind<-df[,11]
       wind[wind<0.1]<-0.1
      # P deposition data
      fid<-open.ncdf(paste("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc",sep=''))
       lonpdep<-as.numeric(get.var.ncdf(fid,'lon'))
       latpdep<-as.numeric(get.var.ncdf(fid,'lat'))
       latx<-which.min(abs(latpdep-latsite))
       lonx<-which.min(abs(lonpdep-lonsite))
       p_dummy<-get.var.ncdf(fid,'pdep')[lonx,latx]/365
       pdep<-numeric(length(swdown))
       pdep[]<-p_dummy
      close.ncdf(fid)
      # N deposition data
    year_out=co2=co2_dC13=co2_DC14=noy=nhx<-numeric(length(swdown))
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      if(year<2017) {
         fid<-open.ncdf(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc",sep=''))
      } else {
         fid<-open.ncdf(paste("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_2016.nc",sep=''))
      }
       lonndep<-get.var.ncdf(fid,'lon')
       latndep<-get.var.ncdf(fid,'lat')
       latx<-which.min(abs(latndep-as.numeric(latsite)))
       lonx<-which.min(abs(lonndep-as.numeric(lonsite)))
       nhx_dummy<-get.var.ncdf(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
       noy_dummy<-get.var.ncdf(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
       nhx[1:(dpm[1]*48)+365*48*(year-startyear)]<-rep(nhx_dummy[1],dpm[1]*48)
       noy[1:(dpm[1]*48)+365*48*(year-startyear)]<-rep(noy_dummy[1],dpm[1]*48)
       for (m in 2:12) {
         nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(nhx_dummy[m],dpm[m]*48)
         noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(noy_dummy[m],dpm[m]*48)
       }
      close.ncdf(fid)
      #CO2
      if(year<2017) {
        co2[1:17520+365*48*(year-startyear)]<-rep(co2y[co2y[,1]==year,2],17520)
      } else {
        co2[1:17520+365*48*(year-startyear)]<-rep(co2y[co2y[,1]==2016,2],17520)+2*(year-2016)
      }
      if(year<2016) {
        co2_dC13[1:17520+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],17520)
      } else {
        co2_dC13[1:17520+365*48*(year-startyear)]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==2015,2],17520)-0.02*(year-2015)
      }
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      if(year<2016) {
        co2_DC14[1:17520+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],17520)
      } else {
        co2_DC14[1:17520+365*48*(year-startyear)]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==2015,C14lat],17520)-5*(year-2015)
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
     #generate site list
#     pftid<-which.max(allsiteinfo[allsitename==site_list[s],15:27])
#     if(pftid==2) pft="TrBE"
#     if(pftid==3) pft="TrBR"
#     if(pftid==4) pft="TeNE"
#     if(pftid==5) pft="TeBE"
#     if(pftid==6) pft="TeBS"
#     if(pftid==7) pft="BNE"
#     if(pftid==8) pft="BBS"
#     if(pftid==9) pft="BNS"
#     if(pftid==10) pft="TeH"
#     if(pftid==11) pft="TrH"
#     if(pftid==12) pft="TeC"
#     if(pftid==13) pft="TrC"
     #slm_site_list[s,]<-unlist(c(site_list[s],lonsite,latsite,pft,startyear,endyear,allsiteinfo[allsitename==site_list[s],4:14]))
  #} else {
  #     print('site not available')
  #}
 
}

#colnames(slm_site_list)<-header_site_list
#write.table(slm_site_list,file=paste(opath,'all_site_list.dat',sep=''),quote=FALSE,row.names=FALSE)
