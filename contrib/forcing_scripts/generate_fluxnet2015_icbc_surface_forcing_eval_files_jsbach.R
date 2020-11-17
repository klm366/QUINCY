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

##
## creates forcing (and surface [ic/bc] + evaluation) files for FLUXNET site-level runs with JSBACH
##
## original script by Juergen Knauer creating files for JSBACH3: paste0(workDirMain, "data_by_jknauer/create_surface_forcing_files_fluxnet.r")
## 
## this version is tailored for producing JSBACH met-forcing files that can be used for creating QUINCY met-forcing files
## the creation of surface and evaluation files is not tested and maybe needs some work, see e.g. cover_type definition
##
## notes:
##   - the surface file creation works for the three sites set as default, but values for e.g. soil depth I do not have for other sites
##   - hence, skipping the surface (ic/bc) file creation, as these files are not needed for QUINCY anyway
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_fluxnet2015_icbc_surface_forcing_eval_files_jsbach.R
##
## met-forcing is created for all sites listed in ./generate_fluxnet2015_site_list.dat by default (see siteListFile)
##

############################## info from Juergens Script docu:
##### Open questions:
## is max_moist calculation correct?
## why do calculated values agree so badly with the values in the Email from Sönke?
##### TODO
## enable case of 2 or more tiles
###################################################################################################

## get the dir this script is located
args      <- commandArgs(trailingOnly = F)  
scriptDir <- paste0(normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)]))), "/") 
print(paste0("This script is located in: ", scriptDir))

## libraries
library(ncdf4)
library(bigleaf)   # functions: aerodynamic.conductance() & surface.conductance() & VPD.to.q() & VPD.to.rH()


###################################################################################################

## definition of custom functions used here (implemented by Jürgen) and previously sourced from:
  # source("JSBACH_postprocessing_helper_functions.r")
  # source("site_setup_functions.r")


## test if year is a leap year
leap.year <- function(year){
  ly <- FALSE
  if (year %% 4 == 0 && year %% 100 != 0 || year %% 400 == 0){
    ly <- TRUE
  }
  return(ly)
}

## timestep: halfhourly, hourly, or daily
get.month.day.from.year <- function(startyear, endyear, timestep){
  
  if(timestep == "daily"){
    mf <- 1
  } else if (timestep == "hourly"){
    mf <- 24
  } else if (timestep == "halfhourly"){
    mf <- 48
  } else {
    stop("'timestep' has to be one of 'daily', 'hourly', 'halfhourly'")
  }
  
  mly <- c(31,29,31,30,31,30,31,31,30,31,30,31)   # number of days per month in a leap year
  mny <- c(31,28,31,30,31,30,31,31,30,31,30,31)   ## number of days per month in a normal year
  
  month = day = doy = numeric()
  for (yr in startyear:endyear){
    
    ly <- leap.year(yr)
    
    if (ly){
      month <- c(month,rep(c(1:12),(mly*mf)))
      day   <- c(day,unlist(sapply(mly, function(x) sort(rep(1:x,mf)))))
      doy   <- c(doy,sort(rep(c(1:366),mf)))
    } else {
      month <- c(month,rep(c(1:12),(mny*mf)))
      day   <- c(day,unlist(sapply(mny, function(x) sort(rep(1:x,mf)))))
      doy   <- c(doy,sort(rep(c(1:365),mf)))
    }
  }
  return(data.frame("month"=month,"day"=day,"doy"=doy))
}

## get solar elevation angle in rad & deg
SolarElevationAngle <- function(doy, hour, lat, lon, useSolartime=FALSE){
  
  SolTime <- if(useSolartime){ 
    # Correction for local time and equation of time
    # Hour + LocTime_h + EqTime_h
    stop("not implemented yet!")
  } else {
    #! Note: For reproducing values close to Fluxnet Rg_pot which is without local time and eq of time correction
    #! (CEIP is even different)
    warning('Solar position calculated without correction for local time and equation of time.')
    hour
  }

  # Conversion to radians
  SolTime_rad <- (SolTime - 12) * pi / 12.0
  # Correction for solar time < -pi to positive, important for SolAzim_rad below
  SolTime_rad <- ifelse(SolTime_rad < -pi, SolTime_rad+2*pi, SolTime_rad)
  
  
  # Fractional year in radians
  FracYear_rad <- 2 * pi * (doy-1) / 365.24  
  
  # Solar declination in radians, accounting for the earth axis tilt
  SolDecl_rad <- (0.33281-22.984*cos(FracYear_rad) - 0.34990*cos(2*FracYear_rad) - 0.13980*cos(3*FracYear_rad)
                    + 3.7872*sin(FracYear_rad) + 0.03205*sin(2*FracYear_rad) + 0.07187*sin(3*FracYear_rad))/180*pi 
                                
   # Solar elevation (vertical, zenithal angle) in radians with zero for horizon
  SolElev_rad <-  asin(sin(SolDecl_rad) * sin(lat/180*pi)
                      + cos(SolDecl_rad) * cos(lat/180*pi) * cos(SolTime_rad))

  SolElev_deg <- SolElev_rad * 180/pi
  
  return(data.frame("SolElev_rad"=SolElev_rad,"SolElev_deg"=SolElev_deg))
}

## get extraterrestrial solar radiation
ExtRadiation <- function(Doy){    
  # Calculate extraterrestrial solar radiation after Lanini, 2010 (Master thesis, Bern University)
  # Fractional year in radians
  FracYear_rad.V.n <- 2*pi*(Doy-1) /365.24
  
  # Total solar irradiance
  SolarIrr_Wm2.c <- 1366.1 #W/m-2
  
  #Eccentricity correction
  ExtRadiation.V.n <- SolarIrr_Wm2.c * (1.00011 + 0.034221*cos(FracYear_rad.V.n) + 0.00128*sin(FracYear_rad.V.n)
                                        + 0.000719*cos(2*FracYear_rad.V.n) + 0.000077*sin(2*FracYear_rad.V.n))
  attr(ExtRadiation.V.n, 'varnames') <- 'ExtRad'
  attr(ExtRadiation.V.n, 'units') <- 'W_m-2'
  
  return(ExtRadiation.V.n)
  
}

##  get potential radiation
PotRadiation <- function(Doy,Hour,Lat,Long,SolElev_deg){
  ##use Solartime (where noon is exactly at 12:00). Set this to FALSE to directly use local winter time
  ## original function from: https://github.com/bgctw/REddyProc/blob/master/R/GeoFunctions.R (27.07.2017)
  
  # Calculate potential radiation from solar elevation and extraterrestrial solar radiation
  SolElev_rad  <- SolElev_deg *pi/180
  ExtRadiation <- ExtRadiation(Doy)
  PotRadiation <- ifelse(SolElev_rad <= 0, 0, ExtRadiation * sin(SolElev_rad) )
  
  return(PotRadiation)
}

## get soil parameters from texture
funct_soilParamFromTexture <- function(clay_fraction, sand_fraction, silt_fraction, rootdepth){
  
  if(sum(clay_fraction,sand_fraction,silt_fraction) != 1){
    stop("function funct_soilParamFromTexture(): grain size distribution has to sum up to 1.0!")
  }
  
  clay <- clay_fraction * 100
  sand <- sand_fraction * 100
  silt <- silt_fraction * 100
  
  #Saxton et al. 1986, Table 2
  a=-4.396;b=-0.0715;c=-4.88*1e-4;d=-4.285*1e-5;e=-3.14;f=-2.22*1e-3;
  g=-3.484*1e-5;h=0.332;m=-0.108;n=0.341;j=-7.251*1e-4;k=0.1276;p=12.012;
  q=-7.55*1e-2;r=-3.895;t=3.671*1e-2;u=-0.1103;v=8.7546*1e-4;
  
  A=exp(a+b*clay+c*(sand^2)+d*(sand^2)*clay)*100.;
  B=e+f*(clay^2)+g*(sand^2)*clay;
  
  ## field capacity = 60 hPa -> 6 kPa
  vfc=exp(log(6/A)/B);
  
  ## permanent wilting point = 10^4.2 ha
  pot_pwp=(10^4.2)/10;
  vpwp=exp(log(pot_pwp/A)/B);
  
  ## available water content (more correctly field capacity of the entire soil column)
  awc_m=vfc
  awc=rootdepth*vfc
  
  ## total pore space
  vpor = h + j*sand+k*log10(clay)
  
  ## mpot
  mpot <- 100*(m+n*vpor)/1000
  
  ## hydraulic conductivity
  ksat = 2.778 * 1e-6 * exp(p + q*sand + (r + t*sand + u*clay + v*clay^2) / vpor)
  
  ## Clapp Hornberger 1978, Table2
  t1 <- c(0.03,0.06,0.09,0.14,0.19,0.28,0.34,0.34,0.43,0.49,0.63)
  t2 <- c(4.05,4.38,4.9,5.3,5.39,7.12,7.75,8.52,10.4,10.4,11.4)
  i=1
  clapp_b <- t2[i]
  while (i < 10){
    if((t1[i]+t1[i+1])/2>=clay_fraction){
      clapp_b=t2[i]
      break
    }
    i=i+1
  }
  if (clay_fraction >=(t1[10] + t1[11])/2){
    clapp_b=t2[11]
  }
  
  
  ## USDA soil triangle
  stype <- 'm'
  if(clay_fraction>0.35){
    stype='f'
  }
  if (clay_fraction<0.35 & sand_fraction < 0.15){
    stype='mf'
  }
  if (clay_fraction>0.6){
    stype='vf'
  }
  if (clay_fraction<0.2 & sand_fraction > 0.65 ){
    stype='c'
  }
  
  #SSSA definition
  stype='L'
  if (clay_fraction>0.4){
    stype='C'
    if (sand_fraction>0.8){ 
      stype='sC'
    }
    if (silt_fraction>0.4){ 
      stype='siC'
    }
  } else {
    if (clay_fraction>0.25){
      stype='cL'
      if(sand_fraction<0.5){ 
        stype='sicL'
      } 
      if(sand_fraction>0.7){ 
        stype='scL'
      }   
    } else {
      if (silt_fraction>0.5){
        stype='siL'
      } else if (silt_fraction>0.28 & sand_fraction < 0.7){
        stype='L'
      } else if (sand_fraction>0.9){ 
        stype='S'
      } else if (sand_fraction>0.75){ 
        stype='lS'
      } else { 
        stype='sL'
      }
    } 
  }
  
  ##Van Genuchten M according to Schaap et al. SSSA 1998, 62, 547ff. 
  if (stype == 'S'){
    logn=0.462
  } else if (stype == "lS"){
    logn=0.194
  } else if (stype == 'sL'){
    logn=0.15
  } else if (stype == 'L'){
    logn=0.195
  } else if (stype == 'siL'){
    logn=0.26
  } else if (stype == 'scL'){
    logn=0.132
  } else if (stype == 'cL'){
    logn=0.188
  } else if (stype == 'sicL'){
    logn=0.24
  } else if (stype == 'C'){
    logn=0.107
  }
  
  n    <- exp(logn)
  m    <- 1-1/n
  spsi <- m/(1-m)
  
  return(c("field_capacity"=vfc,
           "max_moist"=awc,
           "wilting_point"=vpwp,
           "soil_porosity"=vpor,
           "moisture_pot"=mpot,
           "hyd_cond_sat"=ksat,
           "bclapp"=clapp_b,
           "pore_size_index"=spsi))
}  # funct_soilParamFromTexture


###################################################################################################

## define main dirs
#debug baseDir     <- "/home/jengel/remote/rt-bgi/"
baseDir     <- "/Net/Groups/BGI/"
workDirMain <- paste0(baseDir,"work_2/QUINCY/SLM/Data/fluxnet2015_jsbach_input/")


## forcing data input dir: FLUXNET2015 Berkeley dataset
# 1/2 hourly
fluxnetForcingHHDir   <- paste0(baseDir, "data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/Data/HH/2016_11/")
# monthly
fluxnetForcingMMDir   <- paste0(baseDir, "data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/Data/MM/2016_11/")

## forcing info - what sites do have only hourly data
#  for some sites only HR (hourly) data is available (altough they are saved in the HH dir) - ofr these Uli created a HH time series anyway by leaving every timestep to "missing data"
#  this dir contains two .csv files per site
fluxForcHRListDir     <- paste0(baseDir, "data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/orgdata/org2016_11/csv/HR/")

## output dir
surfaceOutDir         <- paste0(workDirMain, "ic_bc_surface/")
                          # original dir: model/MPIESM/input/jsbach/FLUXNET/surface/v2/
modelForcingOutDir    <- paste0(workDirMain, "met_forcing/")
                          # original dir: model/MPIESM/input/jsbach/FLUXNET/halfhourly/FLUXNET2015/
modelEvalOutDir       <- paste0(workDirMain, "evaluation/")
                          # original dir: model/MPIESM/evaluation/jsbach/FLUXNET/halfhourly/FLUXNET2015/

## template surface (output) file
template_surface      <- paste0(workDirMain, "data_by_jknauer/AT-Neu.1_tiles.5_layers.surface_template.nc")
                          # original file from (2019-11-20): people/jknauer/mesophyll_conductance/JSBACH/new_sites_setup/AT-Neu.1_tiles.5_layers.surface_template.nc

d13C_data             <- read.table(paste0(workDirMain, "data_by_jknauer/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt"), header=T, skip=1)
                          # original file from (2019-11-20): people/jknauer/mesophyll_conductance/isotopes/global_delta13C_data/delta13C_in_air_input4MIPs_GM_1850-2015_V1.txt 
CO2_data              <- read.table(paste0(workDirMain, "data_by_jknauer/GCP_co2_global.dat"), header=F)
                          # original file from (2019-11-20): people/jknauer/mesophyll_conductance/JSBACH/new_sites_setup/GCP_co2_global.dat

## site list file with surface properties
#used for ic/bc surface file creation: siteinfo              <- read.table(paste0(workDirMain, "data_by_jknauer/new_fluxnet_sites_surface_properties.txt"), header=T, sep="\t")


## soil layers: valid for all sites
levels     <- c(0.065,0.25,0.91,2.91,5.7)        # soil layer boundaries
depths     <- (c(levels,0) - c(0,levels))[1:5]   # depths of soil layers

## define sites for which forcing should be created
siteListFile <- read.table(paste0(scriptDir, "generate_fluxnet2015_site_list.dat"), header=T, sep="", stringsAsFactors=F)
sites        <- siteListFile$site_id



## cd into working dir
setwd(workDirMain)


###################################################################################################
## I) Surface files a.k.a. ic/bc files (init/boundary conditions)
###################################################################################################

# switch on/off
createSurface <- F   # T/F

### Warning: has only been tested for 1-tile surface files!!!!

if (createSurface){
for (s in 1:length(sites)){
  
  print(paste0("surface files - starting site: ",sites[s]))
  
  cover_type    <- siteinfo[siteinfo$site==sites[s],"cover_type"]    # <- quincy PFT number 1-8
  elevation     <- siteinfo[siteinfo$site==sites[s],"elevation"]
  LAI           <- siteinfo[siteinfo$site==sites[s],"lai"]
  latitude      <- siteinfo[siteinfo$site==sites[s],"lat"]
  longitude     <- siteinfo[siteinfo$site==sites[s],"lon"]
  soil_depth    <- siteinfo[siteinfo$site==sites[s],"soil_depth"]
  root_depth    <- siteinfo[siteinfo$site==sites[s],"root_depth"]
  fraction_clay <- siteinfo[siteinfo$site==sites[s],"fraction_clay"]
  fraction_sand <- siteinfo[siteinfo$site==sites[s],"fraction_sand"]
  fraction_silt <- siteinfo[siteinfo$site==sites[s],"fraction_silt"]
  
  ## if necessary, convert longitude to degrees east
  if (longitude < 0){   # convert lon to degrees east!
    longitude <- 360 + longitude 
  }
  
  ## get mean monthly Tair (which is needed for the surface file...)
  # get from filename the year start and end (e.g. DE-Hai.HH.2000.2012.nc)
  startyear <- substr(list.files(fluxnetForcingHHDir, pattern=sites[s])[1],11,14)
  endyear   <- substr(list.files(fluxnetForcingHHDir, pattern=sites[s])[1],16,19)
  nryears   <- length(startyear:endyear)
  
  forcing_monthly   <- nc_open(paste0(fluxnetForcingMMDir,sites[s],".MM.",startyear,".",endyear,".nc"))
  Tair_monthly      <- cbind(rep(1:12,nryears),ncvar_get(forcing_monthly,"TA_F"))
  Tair_monthly_mean <- round(aggregate(Tair_monthly,by=list(Tair_monthly[,1]),mean,na.rm=T)[,3],3)
  
  
  ## load template surface file and calculate soil parameters and layer moist based on grain size distribution
  templateSurf <- nc_open(template_surface, write=TRUE)
  
  soilpars    <- funct_soilParamFromTexture(fraction_clay, fraction_sand, fraction_silt, root_depth)
  layer_moist <- depths * soilpars["field_capacity"] * 0.5
  
  
  ncvar_put(templateSurf,"bclapp",rep(soilpars["bclapp"],2))
  ncvar_put(templateSurf,"cover_fract",c(1.0,0))
  ncvar_put(templateSurf,"cover_type",c(cover_type,0))
  ncvar_put(templateSurf,"elevation",rep(elevation,2))
  ncvar_put(templateSurf,"forest_fract",rep(1.0,2))
  ncvar_put(templateSurf,"fraction_clay",rep(fraction_clay,2))
  ncvar_put(templateSurf,"fraction_sand",rep(fraction_sand,2))
  ncvar_put(templateSurf,"fraction_silt",rep(fraction_silt,2))
  ncvar_put(templateSurf,"hyd_cond_sat",rep(soilpars["hyd_cond_sat"],2))
  ncvar_put(templateSurf,"init_moist",rep(soilpars["max_moist"],2))              # = max_moist (field capacity of the entire soil column (m))
  ncvar_put(templateSurf,"lai_clim",matrix(LAI,nrow=2,ncol=12))
  ncvar_put(templateSurf,"lat",c(latitude,latitude + 1))
  ncvar_put(templateSurf,"lon",longitude)
  ncvar_put(templateSurf,"layer_moist",rbind(layer_moist,rep(0,5)))
  ncvar_put(templateSurf,"maxmoist",rep(soilpars["max_moist"],2))  
  ncvar_put(templateSurf,"moisture_pot",rep(soilpars["moisture_pot"],2))
  ncvar_put(templateSurf,"pore_size_index",rep(soilpars["pore_size_index"],2))
  ncvar_put(templateSurf,"root_depth",rep(root_depth,2))                         # see email
  ncvar_put(templateSurf,"soil_depth",rep(soil_depth,2))                         # not known, set to 1m
  ncvar_put(templateSurf,"soil_field_cap",rep(soilpars["field_capacity"],2))     # volumetric field capacity (m/m)
  ncvar_put(templateSurf,"soil_porosity",rep(soilpars["soil_porosity"],2))
  ncvar_put(templateSurf,"surf_temp",matrix(c(Tair_monthly_mean,rep(0,12)),nrow=2,ncol=12,byrow=T))  # Tsurf = Tair
  ncvar_put(templateSurf,"wilting_point",rep(soilpars["wilting_point"],2))
  
  ## close file
  nc_close(templateSurf)

  ## copy and rename output surface file
  file.copy(template_surface, paste0(surfaceOutDir,sites[s], ".1_tiles.5_layers.surface.nc"), overwrite=TRUE)

} # end site loop
} # if (createSurface)

##### (comment by Juergen) modify existing surface files manually (maybe it's easier in that case to let it rain more...)
# surf <- nc_open(paste0(surfaceOutDir,"FR-Pue.2_tiles.5_layers.surface.nc"),write=TRUE)




#################################################################################################
## II) Forcing files  (plus eval files, if "create_eval_files == T") 
#################################################################################################

## Settings for forcing (and evaluation) files creation
interpolate          <- T  # if hourly resolution, perform linear interpolation?
create_forcing_files <- T  # create forcing files?
create_eval_files    <- F  # create evaluation files?
calc_Gs_Ga           <- T  # this is used only if create_eval_files = TRUE: calculate Gs and Ga from EC data?


## Important!! 
##  For some sites only HR (hourly) data is available, altough they are saved in the HH (1/2 hourly) dir.
##  For these sites, Uli created a HH time series anyway by leaving every timestep to "missing data".
##  The sites with only HR data are listed in the dir: "$baseDir"data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/orgdata/org2016_11/csv/HR/
## 
##  For this hourly data in halfhourly resolution, it is assumed that every second value is missing, and the first one is missing

## loop over sites and create forcing files (and eval files)
for (s in 1:length(sites)){

  print(" ")
  print("-------------------------")
  print(paste0("creating forcing files - starting site: ", sites[s], " (", s, " of ", length(sites), " sites)"))
  
  ## get info about whether the data is l_HH T/F
  # logical 1/2 data per default == T
  l_HH             <- T
  secondspertstep  <- 3600 / 2
  forcingFileHR    <- list.files(fluxForcHRListDir, pattern=sites[s])
  if (length(forcingFileHR) > 0){
    l_HH <- F
    secondspertstep  <- 3600
    print(" Only hourly data is available for this site")
    if (interpolate){
      print("  Running interpolation HR -> HH")
    } else{
      print("  NOT running interpolation HR -> HH because the Flag interpolate is set to FALSE!")
    }
  }  

  ## get start and end years from file names
  startyear <- substr(list.files(fluxnetForcingHHDir, pattern=sites[s])[1],11,14)
  endyear   <- substr(list.files(fluxnetForcingHHDir, pattern=sites[s])[1],16,19)
  
  cat("forcing file creation - starting site: ",sites[s], fill=TRUE)
  cat("startyear: ",startyear, fill=TRUE)
  cat("endyear: ",endyear, fill=TRUE)

  print(paste0("Opening the original Fluxnet file: ", fluxnetForcingHHDir, sites[s], ".HH.", startyear, ".", endyear, ".nc"))
  fluxForcFile <- nc_open(paste0(fluxnetForcingHHDir, sites[s], ".HH.", startyear, ".", endyear, ".nc"))
  
  ## start forcing file creation
  if (create_forcing_files){
    
    ## lon and lat
    latitude  <- ncvar_get(fluxForcFile, "latitude")
    longitude <- ncvar_get(fluxForcFile, "longitude")
    if (longitude < 0){   # convert lon to degrees east!
      longitude <- 360 + longitude 
      print(paste0("convert lon to degrees east | longitude is now: ", longitude))
    }
    
    ## load additional time information
    year      <- ncvar_get(fluxForcFile,"year")
    month     <- ncvar_get(fluxForcFile,"month")
    day       <- ncvar_get(fluxForcFile,"day")
    hour      <- ncvar_get(fluxForcFile,"hour")
    doy       <- get.month.day.from.year(startyear,endyear,"halfhourly")[,"doy"]
    minute    <- ncvar_get(fluxForcFile,"minute")
    fract_day <- (hour + ((minute - 15) / 60))/24
    
    date      <- strptime(paste(year,month,day,sep="/"),format="%Y/%m/%d")
    gregtime  <- as.numeric(gsub("-","",as.character(date))) + fract_day
    
    ## load data and correct units
    LW  <- ncvar_get(fluxForcFile,"LW_IN_F")
    PA  <- ncvar_get(fluxForcFile,"PA_F") * 1000
    TA  <- ncvar_get(fluxForcFile,"TA_F")
    SW  <- ncvar_get(fluxForcFile,"SW_IN_F")
    WS  <- ncvar_get(fluxForcFile,"WS_F")
    PR  <- ncvar_get(fluxForcFile,"P_F") / secondspertstep
    VPD <- ncvar_get(fluxForcFile,"VPD_F")
    
    LW_qc <- ncvar_get(fluxForcFile,"LW_IN_F_QC")
    PA_qc <- ncvar_get(fluxForcFile,"PA_F_QC")
    TA_qc <- ncvar_get(fluxForcFile,"TA_F_QC")
    SW_qc <- ncvar_get(fluxForcFile,"SW_IN_F_QC")
    WS_qc <- ncvar_get(fluxForcFile,"WS_F_QC")
    PR_qc <- ncvar_get(fluxForcFile,"P_F_QC")
    VPD_qc<- ncvar_get(fluxForcFile,"VPD_F_QC")
    
    ## calculate mpot  
    solar_elev <- SolarElevationAngle(doy=doy,hour=hour,lat=latitude,lon=longitude,
                                      useSolartime=FALSE)[,"SolElev_deg"]
    Rg_pot <- PotRadiation(Doy=doy,Hour=hour,Lat=latitude,Long=longitude,
                           SolElev_deg=solar_elev)  # W m-2
    # if not HH data
    if(!l_HH){
      Rg_pot[seq(1,length(Rg_pot),2)] <- NA
    }
    Rg_pot_qc <- rep(2,length(Rg_pot))

    ## calculate rH
    rH         <- VPD.to.rH(VPD/10,TA)
    rH[rH < 0] <- 0.002
    rH_qc      <- VPD_qc
    qair_qc    <- VPD_qc
    
    ## fill in CO2 and d13C
    d13C <- rep(d13C_data[d13C_data[,1] %in% (year + 0.5),2],table(year))
    CO2  <- rep(CO2_data[CO2_data[,1] %in% year,2],table(year))
  
    ## combine to matrix and add data quality columns
    forcing_vals <- data.frame(cbind(gregtime,VPD,LW,PA,TA,SW,WS,PR,Rg_pot,
                                     rH,d13C,CO2,LW_qc,TA_qc,SW_qc,WS_qc,
                                     PR_qc,VPD_qc,rH_qc,Rg_pot_qc,qair_qc,PA_qc))

    ## cut forcing data frame to good years
    #forcing_vals <- forcing_vals[year %in% c(good_years_start:good_years_end),]

    ## if interpolate == T and no 1/2 hourly data
    pr_c <- which(colnames(forcing_vals) == "PR")
    if (interpolate & !l_HH){
      cat("data are at hourly resolution -> linearly interpolating to half hours!",fill=TRUE)
      missing <- seq(1,nrow(forcing_vals),2)  # assuming that every second HH value is missing, and the first HH one is missing
      
      ## go over all columns (variables) except PR (with column number pr_c)
      for (c in c(2:ncol(forcing_vals))[-(pr_c-1)]){  # assuming that first column is "gregtime"
        forcing_vals[missing,c] <- sapply(missing, function(x) mean(c(forcing_vals[x+1,c],forcing_vals[x-1,c]),na.rm=T))
      } 
      ## special treatment for precipitation (PR)
      # get the HR timesteps 
      precip_events <- which(forcing_vals[,pr_c] > 0)
      forcing_vals[precip_events-1,pr_c] <- forcing_vals[precip_events,pr_c]
      #forcing_vals[precip_events,pr_c]   <- 0.5 * forcing_vals[precip_events,pr_c]
      forcing_vals[is.na(forcing_vals[,pr_c]),pr_c] <- 0
    }
    
    ## increase precpitation for sites that are too dry (for various reasons)
    if (sites[s] == "US-Ne1"){  # increase precip for irrigated sites
      forcing_vals[,pr_c] <- 1.2*forcing_vals[,pr_c]
      print(" ! US-Ne1 - increase precpitation by factor 1.2 as this site is too dry (for various reasons)")
    } else if (sites[s] == "DE-Geb"){
      forcing_vals[,pr_c] <- 1.7*forcing_vals[,pr_c]
      print(" ! DE-Geb - increase precpitation by factor 1.7 as this site is too dry (for various reasons)")
    } else if (sites[s] == "AT-Neu"){
      forcing_vals[,pr_c] <- 2.4*forcing_vals[,pr_c]
      print(" ! AT-Neu - increase precpitation by factor 2.4 as this site is too dry (for various reasons)")
    } else if (sites[s] == "GF-Guy"){
      forcing_vals[,pr_c] <- 1.7*forcing_vals[,pr_c]
      print(" ! GF-Guy - increase precpitation by factor 1.7 as this site is too dry (for various reasons)")
    } 

    ## delete every second value if interpolate is set to FALSE and no 1/2 hourly data
    if (!interpolate & !l_HH){
      print(" !! deleting every second value because interpolate is set to FALSE and data are not HH !!")
      missing <- seq(1,nrow(forcing_vals),2) # assuming that every second value is missing, and the first one is missing
      forcing_vals <- forcing_vals[-missing,]
    }
    
    
    ## fill variables if only a few entries are missing (that is required to make sure JSBACH is correctly working)
    # not the precipitation (PR) column
    interpolate_columns <- c(2:ncol(forcing_vals))[-which(c(2:ncol(forcing_vals)) == pr_c)]

    print("fill variables if only a few entries are missing (excluding PR)")
    cat("in total",nrow(forcing_vals),"observations",fill=T)    
    for (c in interpolate_columns){ 
      missing <- which(is.na(forcing_vals[,c]))
      if (length(missing) > 0){
        cat("filling",length(missing),"values in variable",colnames(forcing_vals)[c],fill=TRUE)
        w <- 2
        still_missing <- TRUE
        while (still_missing){
          missing <- which(is.na(forcing_vals[,c]))
          for (m in missing){
            forcing_vals[m,c] <- mean(forcing_vals[max(c((m-w)),0):min((m+w),nrow(forcing_vals)),c],na.rm=T)
          }
          still_missing <- sum(is.na(forcing_vals[,c])) > 0
          w <- w + 2
        }
      }
    }
  
    
    ## calculate qair from VPD (requires pressure, that's why it comes here)
    forcing_vals$qair <- VPD.to.q(forcing_vals$VPD/10,forcing_vals$TA,forcing_vals$PA/1000)
    forcing_vals$qair[forcing_vals$qair < 0] <- 1e-05
    
    
    ## create ncdf entries (dimensions + variables)
    lat  <- ncdim_def("lat",units="degrees_north",vals=c(latitude,latitude+1),longname="latitude")
    lon  <- ncdim_def("lon",units="degrees_east",vals=longitude,longname="longitude")
    time <- ncdim_def("time",units="day as %Y%m%d.%f",vals=forcing_vals$gregtime,unlim=TRUE,longname="time axis",calendar="proleptic_gregorian")
  
    air_pressure      <- ncvar_def("air_pressure",units="Pa",longname="Air Pressure",dim=list(lon,lat,time),missval=-9999.0)
    air_pressure_flag <- ncvar_def("air_pressure_flag",units="-",longname="Air Pressure Flag",dim=list(lon,lat,time),missval=-9999.0)
    air_temp          <- ncvar_def("air_temp",units="degC",longname="Air Temperature",dim=list(lon,lat,time),missval=-9999.0)
    air_temp_flag     <- ncvar_def("air_temp_flag",units="-",longname="Air Temperature Flag",dim=list(lon,lat,time),missval=-9999.0)
    longwave          <- ncvar_def("longwave",units="W/m2",longname="Longwave downward radiation",dim=list(lon,lat,time),missval=-9999.0)
    longwave_flag     <- ncvar_def("longwave_flag",units="-",longname="Longwave downward radiation Flag",dim=list(lon,lat,time),missval=-9999.0)
    shortwave         <- ncvar_def("shortwave",units="W/m2",longname="Shortwave downward radiation",dim=list(lon,lat,time),missval=-9999.0)
    shortwave_flag    <- ncvar_def("shortwave_flag",units="-",longname="Shortwave downward radiation Flag",dim=list(lon,lat,time),missval=-9999.0)
    mpot              <- ncvar_def("mpot",units="W/m2",longname="Potential shortwave downward radiation",dim=list(lon,lat,time),missval=-9999.0)
    mpot_flag         <- ncvar_def("mpot_flag",units="-",longname="Potential shortwave downward radiation Flag",dim=list(lon,lat,time),missval=-9999.0)
    wspeed            <- ncvar_def("wspeed",units="m/s",longname="Wind speed",dim=list(lon,lat,time),missval=-9999.0)
    wspeed_flag       <- ncvar_def("wspeed_flag",units="-",longname="Wind speed Flag",dim=list(lon,lat,time),missval=-9999.0)
    rel_humidity      <- ncvar_def("rel_humidity",units="%",longname="Relative Humidity",dim=list(lon,lat,time),missval=-9999.0)
    rel_humidity_flag <- ncvar_def("rel_humidity_flag",units="-",longname="Relative Humidity Flag",dim=list(lon,lat,time),missval=-9999.0)
    qair              <- ncvar_def("qair",units="kg/kg",longname="Specific Humidity",dim=list(lon,lat,time),missval=-9999.0)
    qair_flag         <- ncvar_def("qair_flag",units="-",longname="Specific Humidity Flag",dim=list(lon,lat,time),missval=-9999.0)
    precip            <- ncvar_def("precip",units="kg/m2/s",longname="Precipitation (Snow + Rain)",dim=list(lon,lat,time),missval=-9999.0)
    precip_flag       <- ncvar_def("precip_flag",units="-",longname="Precipitation (Snow + Rain) Flag",dim=list(lon,lat,time),missval=-9999.0)
    d13C              <- ncvar_def("d13C",units="permil",longname="atmospheric delta13C",dim=list(lon,lat,time),missval=-9999.0)
    CO2               <- ncvar_def("CO2",units="ppm",longname="atmospheric CO2 concentration",dim=list(lon,lat,time),missval=-9999.0)
  
  
    ### create nc file on disk
    filenameForcing     <- paste0(sites[s],".", startyear,"-", endyear,".forcing.halfhourly.nc")
    filenameForcingWork <- paste0(sites[s],".", startyear,"-", endyear,".forcing.halfhourly",Sys.Date(),".nc")
    print(paste0("Creating new forcing file in work dir: ", filenameForcingWork))
    forcing <- nc_create(filename=filenameForcingWork,
                         vars=list(air_pressure,air_pressure_flag,air_temp,air_temp_flag,longwave,longwave_flag,shortwave,
                                   shortwave_flag,mpot,mpot_flag,wspeed,wspeed_flag,rel_humidity,rel_humidity_flag,
                                   qair,qair_flag,precip,precip_flag,d13C,CO2))
  
  
    ### fill nc file with values (if the value equals the missval, it is automatically set to NA,
    ### so it should be the same to just leave it empty)
    ncvar_put(forcing,air_pressure,vals=forcing_vals[,"PA"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,air_pressure_flag,vals=forcing_vals[,"PA_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,air_temp,vals=forcing_vals[,"TA"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,air_temp_flag,vals=forcing_vals[,"TA_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,longwave,vals=forcing_vals[,"LW"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,longwave_flag,vals=forcing_vals[,"LW_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,shortwave,vals=forcing_vals[,"SW"],start=c(1,1,1),count=c(1,1,-1)) 
    ncvar_put(forcing,shortwave_flag,vals=forcing_vals[,"SW_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,mpot,vals=forcing_vals[,"Rg_pot"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,mpot_flag,vals=forcing_vals[,"Rg_pot_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,wspeed,vals=forcing_vals[,"WS"],start=c(1,1,1),count=c(1,1,-1)) 
    ncvar_put(forcing,wspeed_flag,vals=forcing_vals[,"WS_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,rel_humidity,vals=forcing_vals[,"rH"],start=c(1,1,1),count=c(1,1,-1)) 
    ncvar_put(forcing,rel_humidity_flag,vals=forcing_vals[,"rH_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,qair,vals=forcing_vals[,"qair"],start=c(1,1,1),count=c(1,1,-1)) 
    ncvar_put(forcing,qair_flag,vals=forcing_vals[,"qair_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,precip,vals=forcing_vals[,"PR"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,precip_flag,vals=forcing_vals[,"PR_qc"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,d13C,vals=forcing_vals[,"d13C"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(forcing,CO2,vals=forcing_vals[,"CO2"],start=c(1,1,1),count=c(1,1,-1))
  
    nc_close(forcing)
  
    ## copy result file to the right folder
    print("Copying the newly created forcing file into the met_forcing dir overwriting any previous version & keeping the version in the work dir for manual deletion")
    print(paste0(filenameForcingWork, " -> ", modelForcingOutDir, filenameForcing))
    file.copy(from=filenameForcingWork, to=paste0(modelForcingOutDir, filenameForcing), overwrite=TRUE)
  
  }  # end: if (create_forcing_files)
  
  
  #################################################################################################
  ## III) Observation files (evaluation files)
  #################################################################################################
  if (create_eval_files){

    ## interim version, just to make it work
    ## TODO: the corresponding flag files must be added as well!!
    
    GPP  <- ncvar_get(fluxForcFile,"GPP_NT_VUT_USTAR50")
    H    <- ncvar_get(fluxForcFile,"H_F_MDS")
    LE   <- ncvar_get(fluxForcFile,"LE_F_MDS")
    NEE  <- ncvar_get(fluxForcFile,"NEE_VUT_USTAR50")
    PPFD <- ncvar_get(fluxForcFile,"PPFD_IN")
    Reco <- ncvar_get(fluxForcFile,"RECO_NT_VUT_USTAR50")
    Rn   <- ncvar_get(fluxForcFile,"NETRAD")
    
    if (calc_Gs_Ga){
      Tair  <- ncvar_get(fluxForcFile,"TA_F_MDS")
      Press <- ncvar_get(fluxForcFile,"PA_F")
      VPD   <- ncvar_get(fluxForcFile,"VPD_F_MDS")/10
      wind  <- ncvar_get(fluxForcFile,"WS_F")
      ustar <- ncvar_get(fluxForcFile,"USTAR")
    }

    if (!create_forcing_files){
      ## lon and lat
      latitude  <- ncvar_get(fluxForcFile, "latitude")
      longitude <- ncvar_get(fluxForcFile, "longitude")
      if (longitude < 0){   # convert lon to degrees east!
        longitude <- 360 + longitude 
      }
      
      # load time information
      year  <- ncvar_get(fluxForcFile,"year")
      month <- ncvar_get(fluxForcFile,"month")
      day   <- ncvar_get(fluxForcFile,"day")
      hour  <- ncvar_get(fluxForcFile,"hour")
      doy   <- get.month.day.from.year(startyear,endyear,"halfhourly")[,"doy"]
      minute <- ncvar_get(fluxForcFile,"minute")
      fract_day <- (hour + ((minute - 15) / 60))/24
      
      date <- strptime(paste(year,month,day,sep="/"),format="%Y/%m/%d")
      gregtime <- as.numeric(gsub("-","",as.character(date))) + fract_day
    
    } 
      gregtime <- gregtime[year %in% c(good_years_start:good_years_end)]
      
      
    ## create ncdf entries (dimensions + variables)
    lat  <- ncdim_def("lat",units="degrees_north",vals=latitude,longname="latitude")
    lon  <- ncdim_def("lon",units="degrees_east",vals=longitude,longname="longitude")
    time <- ncdim_def("time",units="day as %Y%m%d.%f",vals=gregtime,unlim=TRUE,longname="time axis",calendar="proleptic_gregorian")
      


    if (!calc_Gs_Ga){
      obs_vals <- data.frame(cbind(year,GPP,H,LE,NEE,PPFD,Reco,Rn))
    } else {
      obs_vals <- data.frame(cbind(year,GPP,H,LE,NEE,PPFD,Reco,Rn,Tair,Press,VPD,wind,ustar))
    }
    
    ## cut evaluation data frame to good years
    obs_vals <- obs_vals[obs_vals[,"year"] %in% c(good_years_start:good_years_end),]
    
        
    ### interpolate, if necessary
    if (interpolate & !l_HH){
      cat("data are at hourly resolution -> linearly interpolating to half hours!",fill=TRUE)
      missing <- seq(1,nrow(obs_vals),2)  # assuming that every second value is missing, and the first one is missing
      
      for (c in c(2:ncol(obs_vals))){  # assuming that first column is "year"
        obs_vals[missing,c] <- sapply(missing, function(x) mean(c(obs_vals[x+1,c],obs_vals[x-1,c]),na.rm=T))
      } 
    }
    
    
    missing_vals <- rep(-9999.0,nrow(obs_vals))
    
    
    ## calculate Gs and Ga, if required
    if (calc_Gs_Ga){
      Ga <- aerodynamic.conductance(obs_vals,pressure="Press",wind="wind",ustar="ustar")[,"Ga_h"]
      Gs <- surface.conductance(obs_vals,pressure="Press",Ga=Ga)[,"Gs_ms"]
    } else {
      Ga <- missing_vals
      Gs <- missing_vals
    }
      
    #### create variable netcdf entries
    fapar_cyclopes      <- ncvar_def("fapar_cyclopes",units="-",longname="fapar_cyclopes",dim=list(lon,lat,time),missval=-9999.0)
    fapar_cyclopes_flag <- ncvar_def("fapar_cyclopes_flag",units="-",longname="fapar_cyclopes_flag",dim=list(lon,lat,time),missval=-9999.0)
    fapar_modis         <- ncvar_def("fapar_modis",units="-",longname="fapar_modis",dim=list(lon,lat,time),missval=-9999.0)
    fapar_modis_flag    <- ncvar_def("fapar_modis_flag",units="-",longname="fapar_modis_flag",dim=list(lon,lat,time),missval=-9999.0)
    fapar_seawifs       <- ncvar_def("fapar_seawifs",units="-",longname="fapar_seawifs",dim=list(lon,lat,time),missval=-9999.0)
    fapar_seawifs_flag  <- ncvar_def("fapar_seawifs_flag",units="-",longname="fapar_seawifs_flag",dim=list(lon,lat,time),missval=-9999.0)
    Gc                  <- ncvar_def("Gc",units="m/s",longname="canopy conductance",dim=list(lon,lat,time),missval=-9999.0)
    Gc_flag             <- ncvar_def("Gc_flag",units="m/s",longname="canopy conductance flag",dim=list(lon,lat,time),missval=-9999.0)
    gsurf               <- ncvar_def("gsurf",units="mmol m-2 s-1",longname="gapfilled canopy conductance",dim=list(lon,lat,time),missval=-9999.0)
    gsurf_flag          <- ncvar_def("gsurf_flag",units="mmol m-2 s-1",longname="gapfilled canopy conductance flag",dim=list(lon,lat,time),missval=-9999.0)
    GPP                 <- ncvar_def("GPP",units="umol m-2 s-1",longname="gross primary productivity",dim=list(lon,lat,time),missval=-9999.0)
    GPP_flag            <- ncvar_def("GPP_flag",units="umol m-2 s-1",longname="gross primary productivity flag",dim=list(lon,lat,time),missval=-9999.0)
    H                   <- ncvar_def("H",units="W m-2",longname="gapfilled sensible heat flux",dim=list(lon,lat,time),missval=-9999.0)
    H_flag              <- ncvar_def("H_flag",units="W m-2",longname="gapfilled sensible heat flux flag",dim=list(lon,lat,time),missval=-9999.0)
    LE                  <- ncvar_def("LE",units="W m-2",longname="gapfilled latent heat flux",dim=list(lon,lat,time),missval=-9999.0)
    LE_flag             <- ncvar_def("LE_flag",units="W m-2",longname="gapfilled latent heat flux flag",dim=list(lon,lat,time),missval=-9999.0)
    NEE                 <- ncvar_def("NEE",units="umol m-2 s-1",longname="gapfilled net ecosystem exchange",dim=list(lon,lat,time),missval=-9999.0)
    NEE_flag            <- ncvar_def("NEE_flag",units="umol m-2 s-1",longname="gapfilled net ecosystem exchange flag",dim=list(lon,lat,time),missval=-9999.0)
    PPFD                <- ncvar_def("PPFD",units="umol m-2 s-1",longname="photosynthetic photon flux density",dim=list(lon,lat,time),missval=-9999.0)
    PPFD_flag           <- ncvar_def("PPFD_flag",units="umol m-2 s-1",longname="photosynthetic photon flux density flag",dim=list(lon,lat,time),missval=-9999.0)
    Reco                <- ncvar_def("Reco",units="umol m-2 s-1",longname="gapfilled ecosystem respiration",dim=list(lon,lat,time),missval=-9999.0)
    Reco_flag           <- ncvar_def("Reco_flag",units="umol m-2 s-1",longname="gapfilled ecosystem respiration flag",dim=list(lon,lat,time),missval=-9999.0)
    Rn                  <- ncvar_def("Rn",units="W m-2",longname="net energy balance",dim=list(lon,lat,time),missval=-9999.0)
    Rn_flag             <- ncvar_def("Rn_flag",units="W m-2",longname="net energy balance flag",dim=list(lon,lat,time),missval=-9999.0)
    
    filename_eval <- paste0(sites[s],".",good_years_start,"-",good_years_end,".obs.halfhourly.nc")
    evaluation <- nc_create(filename=filename_eval,
                            vars=list(fapar_cyclopes,fapar_cyclopes_flag,fapar_modis,fapar_modis_flag,
                                      fapar_seawifs,fapar_seawifs_flag,Gc,Gc_flag,gsurf,gsurf_flag,GPP,
                                      GPP_flag,H,H_flag,LE,LE_flag,NEE,NEE_flag,PPFD,PPFD_flag,Reco,
                                      Reco_flag,Rn,Rn_flag))
    

    ncvar_put(evaluation,fapar_cyclopes,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,fapar_cyclopes_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,fapar_modis,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,fapar_modis_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,fapar_seawifs,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,fapar_seawifs_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Gc,vals=Gs,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Gc_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,gsurf,vals=Gs,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,gsurf_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,GPP,vals=obs_vals[,"GPP"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,GPP_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,H,vals=obs_vals[,"H"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,H_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,LE,vals=obs_vals[,"LE"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,LE_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,NEE,vals=obs_vals[,"NEE"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,NEE_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,PPFD,vals=obs_vals[,"PPFD"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,PPFD_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Reco,vals=obs_vals[,"Reco"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Reco_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Rn,vals=obs_vals[,"Rn"],start=c(1,1,1),count=c(1,1,-1))
    ncvar_put(evaluation,Rn_flag,vals=missing_vals,start=c(1,1,1),count=c(1,1,-1))
    
    nc_close(evaluation)
    
    ## copy result file to the right folder
    file.copy(from=filename_eval,to=paste0(modelEvalOutDir,filename_eval),
              overwrite=TRUE)
  
  } # end: if (create_eval_files)
  
  ## close the Fluxnet2015 Berkely data file (needed for forcing and evaluation file creation)
  nc_close(fluxForcFile)
  
} # end site loop


print(" #####  fin  ##### ")
