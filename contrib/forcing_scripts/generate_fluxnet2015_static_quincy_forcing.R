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
## generate > static < met forcing for Fluxnet sites using the Fluxnet2015 Berkeley dataset
## note! 
##    the input files used here were created by the script: generate_fluxnet2015_icbc_surface_forcing_eval_files_jsbach.R
##    * these files contain two points (2x lat) with data for only point one
##    * these files contain leap years ! which are 1992 1996 2000 2004 ... each 4 years with the February 29th added (see function testLeapYear())
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_fluxnet2015_static_quincy_forcing.R
## 
## note: allsiteinfo is ./generate_fluxnet2015_site_list.dat 
##       some swdown values in IL-Yat are corrected manually
## 
## output (without leap days) may be generated per default in:
##  ./new_fluxnet2015_static_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)

## seed
set.seed(123)

## set in/out path & create out folder
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/fluxnet2015_jsbach_input/met_forcing/"
opath  <- "./new_fluxnet2015_static_met_forcing/"
dir.create(opath, showWarnings=T, recursive=T)

## input file pattern definitions
infile_ext='.forcing.halfhourly.nc'

## read site info file (produced using OCN ?)
allsiteinfo       <- read.table("generate_fluxnet2015_site_list.dat", sep="", header=T, stringsAsFactor=F)  # created manually using generate_fluxnet_site_list.dat
allsitestartyear  <- allsiteinfo[,2]
allsiteendyear    <- allsiteinfo[,3]
allsitename       <- allsiteinfo[,1]
site_list         <- allsitename

doy <- ceiling(seq(1/48,365,1/48))
hod <- rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)

# define output-file headers
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
             'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

# CO2 info files
co2y      <- read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')

# loop over all sites in the site_list
for (s in 1:length(site_list)) {
#for (s in 1:1) {  # for debugging
  print(site_list[s])
  startyear <- allsitestartyear[s]
  endyear   <- allsiteendyear[s]
  nyear     <- endyear-startyear+1
  
  # test input file for this sites actually exists (file example: AT-Neu.2002-2012.forcing.halfhourly.nc)
  inSiteFile      <- paste0(ipath1, site_list[s], '.', startyear, '-', endyear, infile_ext)
  inSiteFileExist <- file.exists(inSiteFile)
  outSiteFile     <- paste0(opath, site_list[s], "_s_", startyear, "-", endyear, '.dat')  # the "_s_" denotes static met-forcing (compared to transient)

  # if input forcing file exists, start creating forcing
  if(inSiteFileExist) {
    if (file.exists(outSiteFile)) file.remove(outSiteFile) 
    # write header
    write.table(t(cbind(header1,header2)), file=outSiteFile, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=' ') 

    # get meteo data from only  the 1st point (i.e., the 1st latitude)
    fid <- nc_open(inSiteFile)
      lonsite <- ncvar_get(fid,"lon")
      latsite <- ncvar_get(fid,"lat")[1]
      dateInt <- as.integer(ncvar_get(fid,"time"))  # get the dates as integer for removing leap days
      
      swdown  <- ncvar_get(fid,"shortwave")[1,]
      lwdown  <- ncvar_get(fid,"longwave")[1,]
      tair    <- ncvar_get(fid,"air_temp")[1,] +273.15   # K = degC + 273.15
      rain    <- ncvar_get(fid,"precip")[1,] *86400      # precip = rain + snow
      qair    <- ncvar_get(fid,"qair")[1,] *1000
      press   <- ncvar_get(fid,"air_pressure")[1,] /100
      wind    <- ncvar_get(fid,"wspeed")[1,]
      ## forcing input files do not contain snowfall information; QUINCY reads the snow, but does not use it; so, simply set all to zero
      snow    <- rain * 0.0 
      ## correct some swdown values in IL-Yat
      if (site_list[s] == "IL-Yat"){
        print(paste0("manual correction of ", length(which(round(swdown,4) == 275.0768)), " swdown values; set all 275.0768 to 0.0"))
        swdown[which(round(swdown,4) == 275.0768)] <- 0.0
      }
    nc_close(fid)

    # remove leap days from meteo data (all Febr 29th)
    leapPositions <- grep("0229", dateInt)
    if(length(leapPositions) > 0){
      print("removing leap days from meteo data")
      swdown  <- swdown[- leapPositions]
      lwdown  <- lwdown[- leapPositions]
      tair    <- tair[- leapPositions]
      rain    <- rain[- leapPositions]
      qair    <- qair[- leapPositions]
      press   <- press[- leapPositions]
      wind    <- wind[- leapPositions]
      snow    <- snow[- leapPositions]
      # remove date vector; not needed hereafter
      remove(dateInt)
    }

    # P deposition data
    fid <- nc_open("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc")
      lonpdep <- as.numeric(ncvar_get(fid,'lon'))
      latpdep <- as.numeric(ncvar_get(fid,'lat'))
      latx    <- which.min(abs(latpdep-latsite))
      lonx    <- which.min(abs(lonpdep-lonsite))
      p_dummy <- ncvar_get(fid,'pdep')[lonx,latx]/365
      pdep    <- numeric(length(swdown))
      pdep[]  <- p_dummy
    nc_close(fid)
    
    # N deposition data
    # create empty vectors based on the length of the 'swdown' meteo vector
    co2=co2_dC13=co2_DC14=noy=nhx <- numeric(length(swdown))
    for (year in c(startyear:endyear)) {
      print(c(site_list[s],year))
      fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc"))
        lonndep   <- ncvar_get(fid,'lon')
        latndep   <- ncvar_get(fid,'lat')
        latx      <- which.min(abs(latndep-as.numeric(latsite)))
        lonx      <- which.min(abs(lonndep-as.numeric(lonsite)))
        nhx_dummy <- ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
        noy_dummy <- ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
        nhx[1:(dpm[1]*48)+365*48*(year-startyear)] <- rep(nhx_dummy[1],dpm[1]*48)
        noy[1:(dpm[1]*48)+365*48*(year-startyear)] <- rep(noy_dummy[1],dpm[1]*48)
        for (m in 2:12) {
          nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(nhx_dummy[m],dpm[m]*48)
          noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(noy_dummy[m],dpm[m]*48)
        }
      nc_close(fid)
      
      #CO2
      co2[1:17520+365*48*(year-startyear)]      <- rep(co2y[co2y[,1]==year,2],17520)
      co2_dC13[1:17520+365*48*(year-startyear)] <- rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],17520)
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      co2_DC14[1:17520+365*48*(year-startyear)] <- rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],17520)
    }
      
    # create time information for output file
    year_out <- numeric(length(swdown))
    for (year in c(startyear:endyear)) year_out[1:17520+365*48*(year-startyear)] <- rep(year,17520)
    doy_out  <- (rep(doy,nyear))
    hod_out  <- (rep(hod,nyear))
    nstep    <- length(doy_out)

    # combine individual variables into one array for output
    forcing<-cbind(year_out[1:nstep],doy_out,hod_out,swdown[1:nstep],lwdown[1:nstep],tair[1:nstep],qair[1:nstep],
                   press[1:nstep],rain[1:nstep],snow[1:nstep],wind[1:nstep],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                   nhx[1:nstep],noy[1:nstep],pdep[1:nstep])

    # write output
    write.table(format(round(forcing,digits=4),nsmall=1), file=outSiteFile, append=TRUE, quote=FALSE,
                col.names=FALSE, row.names=FALSE, sep=" ")

  } else {
    print(paste0("site file not available: ", inSiteFile))
  }
}


## end ##
