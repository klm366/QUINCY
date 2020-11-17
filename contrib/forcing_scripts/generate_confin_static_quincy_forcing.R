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
## generate > static < met forcing for four finish sites - site-set CONFIN
##
## the original forcing input does contain leap years (1996 2000 2004 2008 2012)
##  in the below script all dates "February 29th" are neglected, i.e., if time = YYYY0229.x
##
## note! data availability differs across sites
##       see object 'siteDataPeriods' for the dates
##       I got these data periods by checking over which periods the variable 'shortwave' is available (i.e., < 5000)
## 
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_confin_static_quincy_forcing.R
## 
## note: * allsiteinfo is ./generate_confin_site_list.dat  
##         (manually created from BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/site_descriptions.txt)
##       * also the actual site info file (confin_sites_1991-2014_list.dat) was created manually using the FI-HYY data from FLUXNET3 site info file
##         (lon/lat from /Net/Groups/BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/site_descriptions.txt)
##       * this site info file needs (manual) modification when a new dataset with more recent data is processed
## 
## output may be generated per default in:
##  ./new_confin_static_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)


## create out folder & set in/out path
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/"
opath  <- "./new_confin_static_met_forcing/"
dir.create(opath , showWarnings=T, recursive=T)

## define & load climate-forcing input file
#  contains the data for all four sites (ut limited to 2003-2013 for swdown & lwdown)
inForcFile <- "climate_4sites_1996-2014_with_leap_years.nc"

## message
print("Generating statci met-forcing for four finish conifers sites (CONFIN)")

## read site info file (manually created from BGI/people/szaehle/Projects/QUINCY/GenerateMANP_Forcing_2018/site_info.dat)
allsiteinfo       <- read.table("generate_confin_site_list.dat", stringsAsFactor=F, header=T)
allsitename       <- allsiteinfo[,1]
allsitestartyear  <- allsiteinfo[,5]
allsiteendyear    <- allsiteinfo[,6]
site_list         <- allsitename

## site index in the met-forcing input file, according to the lon/lat
# lon  23.9597, 24.2430, 24.2948, 26.6378
# lat  60.6418, 67.9872, 61.8474, 67.3619
# site FI-Let,  FI-Ken,  FI-Hyy,  FI-Sod
# site sequence in 'site_list' (i.e., in generate_confin_site_list.dat): "FI-Hyy", "FI-Ken", "FI-Let", "FI-Sod")
siteIndexForcing <- data.frame("FI-Hyy"=3, 
                               "FI-Ken"=2, 
                               "FI-Let"=1, 
                               "FI-Sod"=4)

## meteo data: site specific definition of the period (start:last) where meteo data is available (site sequence identical to 'site_list')
siteDataPeriods <- data.frame("FI-Hyy"=as.numeric(c(17569,  333120)),   # 1997 - 2014
                              "FI-Ken"=as.numeric(c(122737, 315600)),   # 2003 - 2013
                              "FI-Let"=as.numeric(c(227953, 298080)),   # 2009 - 2012 | note, many data points were gap-filled
                              "FI-Sod"=as.numeric(c(87697,  333120)))   # 2001 - 2014

## define time variables
doy=ceiling(seq(1/48,365,1/48))
hod=rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)  # days per month

## number of forcing variables and header names in output file
forcVar <- 17  # number of forcing variables including time (see header1 & 2)
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
           'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

## CO2 input data
co2y      <- read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')


## loop over sites (FI-Hyy, FI-Ken, FI-Let, FI-Sod)
for (s in 1:length(site_list)) {
#for (s in 1:1) {
  
  # message
  print(paste0("start creating met-forcing for site: ", site_list[s]))

  # data availability period is identical for all four sites
  startyear <- allsitestartyear[s]
  endyear   <- allsiteendyear[s]
  nyear     <- endyear-startyear+1
  
  # generated output file
  outSiteFile <- paste0(opath, site_list[s],"_s_",startyear,"-",endyear,'.dat')  # the "_s_" denotes static met-forcing (compared to transient)
  # remove old output file
  if (file.exists(outSiteFile)) file.remove(outSiteFile) 
  
  # write header
  write.table(t(cbind(header1,header2)), file=outSiteFile, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=' ') 

  # define met-forcing input file
  inSiteFile <- paste0(ipath1, inForcFile)

  # get lon/lat from site info file
  lonsite <- allsiteinfo[s,2]
  latsite <- allsiteinfo[s,3]
  
  # get index of this site in the met-forcing file
  sIDForc    <- as.integer(siteIndexForcing[s])

  # get meteo data (data are available for different periods across the four sites)
  dataPeriod <- c(siteDataPeriods[1,s], siteDataPeriods[2,s])  # element range in data vector denoting first and last data point of this period

  fid <- nc_open(inSiteFile)
    dateInt <- as.integer(ncvar_get(fid,"time"))[dataPeriod[1]:dataPeriod[2]]  # get the dates as integer for removing leap days
    
    swdown  <- ncvar_get(fid,"shortwave")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]]
    lwdown  <- ncvar_get(fid,"longwave")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]]
    tair    <- ncvar_get(fid,"air_temp")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]] +273.15   # K = degC + 273.15
    rain    <- ncvar_get(fid,"precip")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]] *86400      # precip = rain (snow not included)
    qair    <- ncvar_get(fid,"qair")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]] *1000
    press   <- ncvar_get(fid,"air_pressure")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]] /100
    wind    <- ncvar_get(fid,"wspeed")[sIDForc, sIDForc, dataPeriod[1]:dataPeriod[2]]
    # some corrections and error catching
    wind[wind < 0.1] <- 0.1   # may not be smaller 0.1 for model stability
    snow    <- rain * 0.0     # forcing input files do not contain snowfall information; QUINCY reads the snow, but does not use it; so, simply set all to zero
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

  # create empty vectors based on the length of the 'swdown' meteo vector
  co2=co2_dC13=co2_DC14=noy=nhx <- numeric(length(swdown))
  # loop over forcing years
  for (year in c(startyear:endyear)) {
    print(c(site_list[s],year))

    ## N deposition
    #  using RCP 8.5 data throughout
    fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_RCP8.5_",year,".nc"))
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
                 
  # write output append yearly data (header created above)
  write.table(format(round(forcing,digits=4),nsmall=1), file=outSiteFile, append=TRUE, quote=FALSE,
              col.names=FALSE, row.names=FALSE, sep=" ")
}


## message
print("Done !")


## end ##
