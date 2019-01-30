#!/usr/bin/env Rscript

#########################################################################################
##
## generate > transient < met forcing for FERT sites using weatherdump v2
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_fert_transient_quincy_forcing.R
## 
## note: allsiteinfo is ./generate_fert_siteinfo_noheader.dat
##       ! only the control sites are created, the high/low N sites may be created by Silvia
## 
## output may be generated per default in:
##  ./new_fert_transient_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)

## create out folder & set in/out path
dir.create("./new_fert_transient_met_forcing", showWarnings=T, recursive=T)
ipath1="/Net/Groups/BGI/people/szaehle/Projects/QUINCY/GenerateCRUForcing_2017/Results/"
ipath2="/Weather/"              # relative path within each site folder
infile_prefix='weather_dump_'
infile_ext='.nc'
opath="./new_fert_transient_met_forcing/"

## forcing period
startyear <- 1901
endyear   <- 2015

## site info
allsiteinfo <- read.table("generate_fert_siteinfo_noheader.dat", stringsAsFactor=F)
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


## loop over all FERT sites and create forcing files
for (s in 1:length(site_list)) {
print(paste0("site: ", site_list[s], "  number ", s, " of ", length(site_list)))
  # check if forcing data are available for this site
  stat3   <- file.exists(paste0(ipath1,site_list[s],ipath2,infile_prefix,1994,infile_ext))
  
  # create name of output file => adding "-CON" to site name for "control"
  outname <- paste0(opath,site_list[s],"-CON_",startyear,"-",endyear,'.dat')
  
  # if forcing data is available
  if(stat3) {
    # remove prev. forcing file 
    if (file.exists(outname)) file.remove(outname) 
    # write header
    write.table(t(cbind(header1,header2)), file=outname, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=' ') 

    # loop over years and write forcing to 'outname'
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
      
      ## create forcing array
      forcing<-cbind(year_out,doy,hod,swdown,lwdown,tair,qair,press,rain,snow,wind,co2,co2_dC13,co2_DC14,nhx,noy,pdep)
      
      ## write forcing to file 
      write.table(format(round(forcing,digits=4),nsmall=1), file=outname,
                  append=TRUE, quote=FALSE, 
                  col.names=FALSE, row.names=FALSE,sep=' ')
    }
  # if NO forcing data is available for this site
  } else {
    print(paste0(" The ", infile_prefix, " file for the site: ", site_list[s], " is not available."))
    # stat2<-unlist(strsplit(stat,split=' '))[11]
    # filestat<-unlist(strsplit(stat2,split="/"))
    # linktarget<-filestat[length(filestat)-1]
    # refname=paste(linktarget,"_",startyear,"-",endyear,'.dat',sep='')
    #    print(c('linking',outname,refname))
    #    system(paste("ln -s",refname,outname))   
  }
}
