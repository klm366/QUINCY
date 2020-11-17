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
## generate > transient < met forcing for four finish sites - site-set CONFIN
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
##  Rscript --vanilla ./generate_confin_transient_quincy_forcing.R
## 
## note: * allsiteinfo is ./generate_confin_site_list.dat  
##         (manually created from BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/site_descriptions.txt)
##       * also the actual site info file (confin_sites_1901-2014_list.dat) was created manually using the FI-HYY data from FLUXNET3 site info file
##         (lon/lat from /Net/Groups/BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/site_descriptions.txt)
##       * this site info file needs (manual) modification when a new dataset with more recent data is processed
## 
## output may be generated per default in:
##  ./new_confin_transient_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)


## seed
set.seed(123)


## create out folder & set in/out path
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/CONFIN_4_finish_sites/"
opath  <- "./new_confin_transient_met_forcing/"
dir.create(opath, showWarnings=T, recursive=T)


## define & load climate-forcing input file
#  contains the data for all four sites (ut limited to 2003-2013 for swdown & lwdown)
inForcFile <- "climate_4sites_1996-2014_with_leap_years.nc"


## message
print("Generating transient met-forcing for four finish conifers sites (CONFIN)")


## read site info file (produced using OCN ?)
allsiteinfo       <- read.table("generate_confin_site_list.dat", stringsAsFactor=F, header=T)
allsitename       <- allsiteinfo[,1]
allsitestartyear  <- allsiteinfo[,5]
allsiteendyear    <- allsiteinfo[,6]
site_list         <- allsitename
allsiteDataYears  <- allsiteendyear - allsitestartyear + 1


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


## create array with random years of the site-forcing for 1901-2008 (108 yrs) per site, to be placed in front of the actual forcing data
## this range of data is not needed for all sites, but the latest year data of confin sites starts (is available) is 2009
periodForcingYears         <- 108
randForcingYears           <- array(NA,dim=c(periodForcingYears, length(site_list)))
colnames(randForcingYears) <- site_list
for (s in 1:length(site_list)){
  randForcingYears[,s] <- sample(1:allsiteDataYears[s], periodForcingYears, replace=T)
}


## create day of year and hour of day and dpm (days per month) default values for output met_forcing files (without leap years)
doy <- ceiling(seq(1/48,365,1/48))
hod <- rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)  


## define output-file headers
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
             'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')


## CO2 info files
co2y      <- read.table('/Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat')
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2015_V1.txt')


## forcing may be created for years 1901 up to 2014 (differs across sites; hence, the endyear is site-specific)
forcing_startyear   <- 1901


## define met-forcing input file
inSiteFile <- paste0(ipath1, inForcFile)


## loop over all sites in the site_list
for (s in 1:length(site_list)) {
#for (s in 1:1) {   # for debugging
  print(site_list[s])
  
  # get the (randomly selected) year of the actual measured forcing for the "transient" forcing period before measured forcing is available
  pickyear <- randForcingYears[,s]

  # start/end year and length of actual measured data
  startyear <- allsitestartyear[s]
  endyear   <- allsiteendyear[s]
  nyear     <- allsiteDataYears[s]
  
  # define output file, e.g., FI-Ken_t_1901-2013.dat  (the '_t_' denotes transient forcing files)
  outname   <- paste0(opath, site_list[s], "_t_", forcing_startyear, "-", endyear, ".dat")

  # replace old output file
  if (file.exists(outname)) file.remove(outname) 

  # write header, space separated
  write.table(t(cbind(header1,header2)), file=outname,
              append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ") 

  # get lon/lat from site info file
  lonsite <- allsiteinfo[s,2]
  latsite <- allsiteinfo[s,3]

  # get index of this site in the met-forcing file
  sIDForc    <- as.integer(siteIndexForcing[s])

  ## get meteo data (data are available for different periods across the four sites)
  dataPeriod <- c(siteDataPeriods[1,s], siteDataPeriods[2,s])  # element range in data vector denoting first and last data point of this period

  # open forcing file and read data
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
  
  ## P deposition data
  fid <- nc_open("/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc")
    lonpdep  <- as.numeric(ncvar_get(fid,'lon'))
    latpdep  <- as.numeric(ncvar_get(fid,'lat'))
    latx     <- which.min(abs(latpdep-latsite))
    lonx     <- which.min(abs(lonpdep-lonsite))
    # P deposition constant value
    p_dummy  <- ncvar_get(fid,'pdep')[lonx,latx]/365
    pdep     <- numeric(length(swdown)) 
    pdep[]   <- p_dummy
  nc_close(fid)
  # create empty vectors used for output data
  yTimeSteps <- 365*24*2              # 30 min time steps in QUINCY
  year_out   <- numeric(yTimeSteps)
  co2        <- numeric(yTimeSteps)
  co2_dC13   <- numeric(yTimeSteps)
  co2_DC14   <- numeric(yTimeSteps)
  noy        <- numeric(yTimeSteps)
  nhx        <- numeric(yTimeSteps)
  ycount     <- 1
  
  # loop over the years from default start of forcing (1901) to the site-specific last years (overall max should be 2014)
  for (year in forcing_startyear:endyear){
    
    ## N deposition
    #  this data covers the period 1850 to 2100
    if (year < 1996){
      # use Ndep data created by Sönke
      fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc"))
    }
    else {
      #  using RCP 8.5 data for all years after 1995
      fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_RCP8.5_",year,".nc"))
    }
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
    
    ## CO2
    co2[1:yTimeSteps]      <- rep(co2y[co2y[,1]==year,2],yTimeSteps)
    co2_dC13[1:yTimeSteps] <- rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],yTimeSteps)
    if(latsite>30) C14lat=2
    if(latsite>-30&latsite<=30) C14lat=3
    if(latsite<=-30) C14lat=4
    co2_DC14[1:yTimeSteps] <- rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],yTimeSteps)
    
    ## year
    year_out[1:yTimeSteps] <- rep(year,yTimeSteps)
    
    ## time information for output as defined above
    doy_out <- doy
    hod_out <- hod
    nstep   <- yTimeSteps
    
    ## get actual year, differentiate between "transient period" and the period with actual fluxnet data
    if(year >= startyear & year <= endyear){
      climyear <- year-startyear+1            # current loop is: for (year in forcing_startyear:endyear)
    } else{
      climyear <- pickyear[ycount]
    }
    
    ## combine individual variables into one array for output (note: nstep == yTimeSteps)
    forcing<-cbind(year_out[1:nstep],doy_out,hod_out,swdown[(1:nstep)+yTimeSteps*(climyear-1)],lwdown[(1:nstep)+yTimeSteps*(climyear-1)],
                   tair[(1:nstep)+yTimeSteps*(climyear-1)],qair[(1:nstep)+yTimeSteps*(climyear-1)],
                   press[(1:nstep)+yTimeSteps*(climyear-1)],rain[(1:nstep)+yTimeSteps*(climyear-1)],
                   snow[(1:nstep)+yTimeSteps*(climyear-1)],wind[(1:nstep)+yTimeSteps*(climyear-1)],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                   nhx[1:nstep],noy[1:nstep],pdep[1:nstep])
    
    ## write output (add another year of data)
    write.table(format(round(forcing,digits=4),nsmall=1), file=outname, append=TRUE, quote=FALSE,
                col.names=FALSE, row.names=FALSE, sep=" ")
    
    ## increase year counter
    ycount <- ycount+1
    ## debug output a the screen
    print(paste0("Created data for site ",site_list[s]," and year ",year," using the fluxnet data from year ",climyear+startyear-1))
  }
} # loop over sites



## end ##
