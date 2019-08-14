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
## generate > transient < met forcing for Fluxnet sites
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_fluxnet_transient_slm_forcing.R
## 
## note: allsiteinfo is ./generate_fluxnet_site_list.dat 
##       some swdown values in IL-Yat are corrected manually
## 
## output may be generated per default in:
##  ./new_fluxnet_transient_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)

## seed
set.seed(123)

## create out folder & set in/out path
dir.create("./new_fluxnet_transient_met_forcing", showWarnings=T, recursive=T)
ipath1 <- "/Net/Groups/BGI/model/OCN/input/point/FLUXNET/v2/"
opath  <- "./new_fluxnet_transient_met_forcing/"
#opath  <- "/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/met_forcing/"

## general extension of the site-specific forcing input file, e.g., DE-Hai.2000-2006.metforcing.nc  in folder ipath1
infile_ext <- ".metforcing.nc"

## read site info file (produced using OCN ?)
allsiteinfo      <- read.table("generate_fluxnet_site_list.dat", stringsAsFactor=F, header=T)  # file is originally from ipath1
allsitename      <- allsiteinfo[,1]
allsitestartyear <- allsiteinfo[,2]
allsiteendyear   <- allsiteinfo[,3]
site_list        <- allsitename
allsiteDataYears <- allsiteendyear - allsitestartyear + 1


## create array with random years of the site-forcing for 1901-2005 (105 yrs) per site
## this range of data is not needed by all sites, but the latest year data of a fluxnet sites starts (is available) is 2006
periodForcingYears         <- 105
randForcingYears           <- array(NA,dim=c(periodForcingYears, length(site_list)))
colnames(randForcingYears) <- site_list
for (s in 1:length(site_list)){
  randForcingYears[,s] <- sample(1:allsiteDataYears[s], periodForcingYears, replace=T)
}


## create day of year and hour of day and dpm (days per month) default values for output met_forcing files
doy <- ceiling(seq(1/48,365,1/48))
hod <- rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)  

## create default headers for output met_forcing files
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
             'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

## get CO2 data
co2y      <- read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

## forcing may be created for years 1901 up to 2006 (differs across sites)
forcing_startyear   <- 1901
forcing_endyear     <- max(allsiteendyear)  # should be 2006


## loop over sites and create per-site met_forcing output
for (s in 1:length(site_list)) {
#for (s in 1:1) {
  print(site_list[s])
  
  pickyear <- randForcingYears[,s]

  # start/end year and length of actual measured data
  startyear <- allsitestartyear[s]
  endyear   <- allsiteendyear[s]
  nyear     <- allsiteDataYears[s]
  
  # check if a site-specific netcdf file exists, e.g., DE-Hai.2000-2006.metforcing.nc
  sDataFile <- paste0(ipath1,site_list[s],".",startyear,"-",endyear,infile_ext)
  stat3     <- file.exists(sDataFile)
  
  # define output file, e.g., DE-Hai_1901-2006.dat or DE-Bay_1901-1998.dat
  outname   <- paste0(opath,site_list[s],"_t_",forcing_startyear,"-",endyear,".dat")

  # if data for this site exists
  if(stat3) {
    # replace old output file
    if (file.exists(outname)) file.remove(outname) 
    # write header, space separated
    write.table(t(cbind(header1,header2)), file=outname,
                append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ") 

    # get met_forcing data
      # meteo data
      fid <- nc_open(sDataFile)
        lonsite <- ncvar_get(fid,"x")
        latsite <- ncvar_get(fid,"y")
        swdown  <- ncvar_get(fid,"SWdown")
        ## correct some swdown values in IL-Yat
        if (site_list[s] == "IL-Yat"){
          print(paste0("manual correction of ", length(which(round(swdown,4) == 275.0768)), " swdown values; set all 275.0768 to 0.0"))
          swdown[which(round(swdown,4) == 275.0768)] <- 0.0
        }
        lwdown  <- ncvar_get(fid,"LWdown")
        tair    <- ncvar_get(fid,"Tair")
        rain    <- ncvar_get(fid,"Rainf")*86400
        snow    <- ncvar_get(fid,"Snowf")*86400
        qair    <- ncvar_get(fid,"Qair")*1000
        press   <- ncvar_get(fid,"PSurf")/100
        wind    <- ncvar_get(fid,"Wind_N")
        wind[wind < 0.1] <- 0.1
      nc_close(fid)
      # P deposition data
      fid <- nc_open(paste("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc",sep=''))
        lonpdep  <- as.numeric(ncvar_get(fid,'lon'))
        latpdep  <- as.numeric(ncvar_get(fid,'lat'))
        latx     <- which.min(abs(latpdep-latsite))
        lonx     <- which.min(abs(lonpdep-lonsite))
        p_dummy  <- ncvar_get(fid,'pdep')[lonx,latx]/365
        pdep     <- numeric(length(swdown))                # P deposition
        pdep[]   <- p_dummy
      nc_close(fid)


    # create empty vectors used for output data
    yTimeSteps <- 365*24*2      # 30 min time steps in SLM
    year_out   <- numeric(yTimeSteps)
    co2        <- numeric(yTimeSteps)
    co2_dC13   <- numeric(yTimeSteps)
    co2_DC14   <- numeric(yTimeSteps)
    noy        <- numeric(yTimeSteps)
    nhx        <- numeric(yTimeSteps)
    ycount     <- 1
    
    # loop over the years from default start of forcing (1901) to the site-specific last years (max 2006)
    for (year in forcing_startyear:endyear){

      # N deposition data
      fid <- nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc"))
        lonndep   <- ncvar_get(fid,'lon')
        latndep   <- ncvar_get(fid,'lat')
        latx      <- which.min(abs(latndep-as.numeric(latsite)))
        lonx      <- which.min(abs(lonndep-as.numeric(lonsite)))
        
        nhx_dummy <- ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
        noy_dummy <- ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
        nhx[1:(dpm[1]*48)] <- rep(nhx_dummy[1],dpm[1]*48)
        noy[1:(dpm[1]*48)] <- rep(noy_dummy[1],dpm[1]*48)
        
        for (m in 2:12) {
          nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(nhx_dummy[m],dpm[m]*48)
          noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48]<-rep(noy_dummy[m],dpm[m]*48)
        }
      nc_close(fid)

      # CO2
      co2[1:yTimeSteps]<-rep(co2y[co2y[,1]==year,2],yTimeSteps)
      co2_dC13[1:yTimeSteps]<-rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],yTimeSteps)
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      co2_DC14[1:yTimeSteps]<-rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],yTimeSteps)

      # year
      year_out[1:yTimeSteps]<-rep(year,yTimeSteps)

      # default output as defined above
      doy_out <- doy
      hod_out <- hod
      nstep   <- yTimeSteps
      
      # get actual year, differentiate between "transient period" and the period with actual fluxnet data
      if(year >= startyear & year <= endyear){
        climyear <- year-startyear+1  # current loop is: for (year in forcing_startyear:endyear)
      } else{
        climyear <- pickyear[ycount]
      }

      # combine individual variables into one array for output
      forcing<-cbind(year_out[1:nstep],doy_out,hod_out,swdown[(1:nstep)+yTimeSteps*(climyear-1)],lwdown[(1:nstep)+yTimeSteps*(climyear-1)],
                     tair[(1:nstep)+yTimeSteps*(climyear-1)],qair[(1:nstep)+yTimeSteps*(climyear-1)],
                     press[(1:nstep)+yTimeSteps*(climyear-1)],rain[(1:nstep)+yTimeSteps*(climyear-1)],
                     snow[(1:nstep)+yTimeSteps*(climyear-1)],wind[(1:nstep)+yTimeSteps*(climyear-1)],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                     nhx[1:nstep],noy[1:nstep],pdep[1:nstep])

      # write output (add another year of data)
      write.table(format(round(forcing,digits=4),nsmall=1), file=outname, append=TRUE, quote=FALSE,
                  col.names=FALSE, row.names=FALSE, sep=" ")

      # increase year counter
      ycount<-ycount+1

      # debug output a the screen
      print(paste0("Created data for site ",site_list[s]," and year ",year," using the fluxnet data from year ",climyear+startyear-1))
    }
  # if no data available for this site do not try to create met_forcing
  } else {
       print('site not available')
  }
} # loop over sites



## end ##
