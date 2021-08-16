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
## generate > transient < met forcing for the EucFace site-set (one site with 12 different metforcing scenarios)
##    forcing is created by combining multiple files differently and using different assumptions about CO2 concentration etc. ...
##    forcing for 12 scenarios is created (number 05-16) in the table at page 4 of:
##      /Net/Groups/BGI/work_2/QUINCY/SLM/Data/EUCFACE/EucFACE_modelling_protocol_2020.pdf
##    met-forcing files will be created with the name EUC_05_.. to EUC_16_... ( see 'scenario_list' below)
##
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_eucface_transient_quincy_forcing.R
## 
## note: allsiteinfo is ./generate_eucface_input_met_data_info.dat created manually by Lin (lyu)
##       it does not contain all 12 resulting met files, but information about the six input met-data files
##       but it is rather used to gather information about the input met-data files !
##       => hence, the variable 'inputDataInfo' is used instead of 'allsiteinfo'
## 
## output may be generated per default in:
##  ./new_eucface_transient_met_forcing/
## 
#########################################################################################


## library
library(ncdf4)

## seed
set.seed(123)

## function to convert VPD to q_air, from:
##   /Net/Groups/BGI/people/szaehle/Projects/FLUXNET/InputPrep/Step2/Scripts/atm_humidity_from_VPD.m
CONVERT_VPD_2_QAIR <- function(Tair,VPD,P){
  #calculates specific humidity from Vapour Pressure Deficit and Air Temp.
  #Input:
  #  Tair Air temperature [°C]
  #  VPD  Vapour Pressure Deficit [Pa]
  #  P Air Pressure [Pa]
  #Output:
  #  Q Specific Humidity [kg/kg]
  #
  #follows Monteith & Unsworth 2008, page 11ff
  #mixing ratio of wet and dry air
  eps=0.622
  # T parameters eq2.27
  Tstar=273.15
  Tstroke=36
  estar=611 # in Pa =.611 kPa
  A=17.27
  
  esat=estar*exp(A*((Tair-Tstar)/(Tair-Tstroke)))
  e = esat - VPD
  e = ifelse (e<0.01*esat,0.01 * esat,e)
  atm_hum = eps * e / ( P - ( 1 -eps ) * e )*1000
  
  return(atm_hum)
}

## set in/out path & create out folder
ipath1 <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/EUCFACE/"
opath  <- "./new_eucface_transient_met_forcing/"
dir.create(opath, showWarnings=T, recursive=T)

## input file pattern definitions
infile_ext='_half_hourly.csv'

## read input met-data info file (produced manually)
inputDataInfo       <- read.table("generate_eucface_input_met_data_info.dat", sep="", header=T, stringsAsFactor=F) 
inputDataName       <- inputDataInfo[,1]
inputDataStartyear  <- inputDataInfo[,2]
inputDataEndyear    <- inputDataInfo[,3]

## site info - it is actually the same site, with 12 different metforcing scenarios
latsite       <- -33.61
lonsite       <- 150.73
startyear     <- 1700
endyear       <- 2069
scenario_list <- c("FIX_AMB_NOP","FIX_ELE_NOP",
                   "FIX_AMB_MDP","FIX_ELE_MDP",
                   "FIX_AMB_HIP","FIX_ELE_HIP",
                   "VAR_AMB_NOP","VAR_ELE_NOP",
                   "VAR_AMB_MDP","VAR_ELE_MDP",
                   "VAR_AMB_HIP","VAR_ELE_HIP")

doy <- ceiling(seq(1/48,365,1/48))
hod <- rep(seq(0,23.5,0.5),365)
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)

# define output-file headers
header1 <- c('year','doy','hour','swvis_srf_down','lw_srf_down','t_air','q_air','press_srf',
             'rain','snow','wind_air','co2_mixing_ratio','co2_dC13','co2DC14','nhx_srf_down','noy_srf_down','p_srf_down')
header2 <- c('-','-','-','Wm-2','Wm-2','K','g/kg','hPa','mm/day','mm/day','m/s','ppm','per-mill','per-mill','mg/m2/day','mg/m2/day','mg/m2/day')

# CO2 info files
# co2y  is included in the original met-data files, see 'ipath1'
co2y_dC13 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/delta13C_in_air_input4MIPs_GM_1850-2019_extrapolated.txt')
co2y_DC14 <- read.table('/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/Atm/Delta14C_in_air_input4MIPs_SHTRNH_1850-2019_extrapolated.txt')

# save the six forcing inputs in a list of DF for further use
#  loop over the files in ipath1 (spinup 1700-1749, observation 1750-2011 & 2012-2019, prediction 2020-2069)
#  based on the names and years in the inputDataInfo file
met_list <- list()
for (m in 1:length(inputDataName)) {
  print(paste0(" reading met-dat input file: ", inputDataName[m], " ", inputDataStartyear[m], "-", inputDataEndyear[m]))
  inputStartyear <- inputDataStartyear[m]
  inputEndyear   <- inputDataEndyear[m]
  if (m==1) {
    filename  <- paste0(inputDataName[m], '_', inputStartyear, '-', inputEndyear)
  } else {
    filename  <- paste0(inputDataName[m], '_', inputStartyear, '_', inputEndyear)
  }
  inSiteFile   <- paste0(ipath1, filename, infile_ext)
  # get input met-data: read separately the header (column names) and the data (fid == file ID, i.e., simply the file)
  fid_colNames <- read.csv(inSiteFile, sep=",", dec=".", header = F, stringsAsFactors = F, nrow = 1)
  fid          <- read.csv(inSiteFile, sep=",", dec=".", header = F, stringsAsFactors = F, skip = 3)
  colnames(fid)       <- unlist(fid_colNames)
  met_list[[m]]       <- fid
  names(met_list)[m]  <- filename
  ## correct the timestep of the forcing input
  #  the forcing input starts each day with 0.5 (i.e., forcing for the 1st half hour)
  #  the default for QUINCY is, that forcing starts each day with 0.0 (also meaning "forcing for the 1st half hour")
  #  as this is just a "naming convention" the 'HOUR' values can simply be reduced by 0.5 w/o any effect on the validity of the forcing for QUINCY  
  met_list[[m]]$HOUR <- met_list[[m]]$HOUR - 0.5
  # cleanup
  remove(fid)
}

# loop over all sites in the scenario_list (12 scenarios)
# and create the forcing input for QUINCY
fid_colNames <- colnames(met_list[[1]])
for (s in 1:length(scenario_list)) {
  #for (s in 1:1) {  # for debugging
  # create scenario file ID
  scenarioID <- s + 4  # first scenario is 05 according to the above mentioned table in EucFACE_modelling_protocol_2020.pdf
  if (scenarioID < 10){
    scenarioFileID <- paste0("EUC_0", scenarioID)
  } else{
    scenarioFileID <- paste0("EUC_",  scenarioID)
  }
  
  # create output filename & path
  outScenarioFile <- paste0(opath, scenarioFileID, "_", scenario_list[s], "_", startyear, "-", endyear, ".dat")
  
  print(paste0(" creating ", outScenarioFile))
  
  # get the correct CO2 value and associate the column name "CO2air" with the column (using the column name from the two files EUC_17*_half_hourly.csv)
  fid0<-met_list[[1]]
  fid1<-met_list[[2]]
  if (grepl("(?=.*FIX)(?=.*ELE)",scenario_list[s],perl=TRUE)) {
    fid2<-met_list[[3]][,-13]
    fid3<-met_list[[5]][,-13]
    colnames(fid2)<-fid_colNames
    colnames(fid3)<-fid_colNames
  } else if (grepl("(?=.*FIX)(?=.*AMB)",scenario_list[s],perl=TRUE)) {
    fid2<-met_list[[3]][,-14]
    fid3<-met_list[[5]][,-14]
    colnames(fid2)<-fid_colNames
    colnames(fid3)<-fid_colNames
  } else if (grepl("(?=.*VAR)(?=.*ELE)",scenario_list[s],perl=TRUE)) {
    fid2<-met_list[[4]][,-13]
    fid3<-met_list[[6]][,-13]
    colnames(fid2)<-fid_colNames
    colnames(fid3)<-fid_colNames
  } else {
    fid2<-met_list[[4]][,-14]
    fid3<-met_list[[6]][,-14]
    colnames(fid2)<-fid_colNames
    colnames(fid3)<-fid_colNames
  }
  
  fid_all <- rbind(fid0,fid1,fid2,fid3)
  rm(fid0,fid1,fid2,fid3)
  
  # if input forcing file exists, start creating forcing
  if(is.data.frame(fid_all)) {
    if (file.exists(outScenarioFile)) file.remove(outScenarioFile) 
    # write header
    write.table(t(cbind(header1,header2)), file=outScenarioFile, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=' ') 
    dateInt   <- as.integer(fid_all[,"DOY"])
    swdown    <- as.numeric(fid_all[,"SWdown"])
    lwdown    <- as.numeric(fid_all[,"LWdown"])
    tair      <- as.numeric(fid_all[,"Tair"])
    rain      <- as.numeric(fid_all[,"Rain"]) *48      # convert from mm/halfhour to mm/day
    vpd       <- as.numeric(fid_all[,"VPD"])
    wind      <- as.numeric(fid_all[,"Wind"])
    wind[wind < 0.1] <- 0.1    # may not be smaller 0.1 for model stability
    ## forcing input files do not contain snowfall information; QUINCY reads the snow, but does not use it; so, simply set all to zero
    snow      <- rain * 0.0 
    pressHLP  <- as.numeric(fid_all[,"PSurf"])
    qair      <- CONVERT_VPD_2_QAIR(tair, vpd, pressHLP)
    press     <- pressHLP / 100
    Ndep      <- as.numeric(fid_all[,"Ndep"])
    co2       <- as.numeric(fid_all[,"CO2air"])
    year_out  <- as.numeric(fid_all[,"YEAR"])
    doy_out   <- as.numeric(fid_all[,"DOY"])
    hod_out   <- as.numeric(fid_all[,"HOUR"])
    rm(fid_all)
    
    # remove leap days from meteo data (all Febr 29th)
    leapPositions <- grep("366", dateInt)
    if(length(leapPositions) > 0){
      print(" removing leap days from meteo data")
      swdown   <- swdown[- leapPositions]
      lwdown   <- lwdown[- leapPositions]
      tair     <- tair[- leapPositions]
      rain     <- rain[- leapPositions]
      vpd      <- vpd[- leapPositions]
      wind     <- wind[- leapPositions]
      snow     <- snow[- leapPositions]
      qair     <- qair[- leapPositions]
      press    <- press[- leapPositions]
      Ndep     <- Ndep[- leapPositions]
      co2      <- co2[- leapPositions]
      year_out <- year_out[- leapPositions]
      doy_out  <- doy_out[- leapPositions]
      hod_out  <- hod_out[- leapPositions]
      # removing 'dateInt' vector (not needed hereafter)
      remove(dateInt)
    }
    nstep <- length(doy_out)
    
    # lon (from the jsbach metforcing files) needs to be converted back to the range [-180,180]
    # it was converted to degrees east in the generate_fluxnet2015_icbc_surface_forcing_eval_files_jsbach.R
    if (lonsite > 180){ 
      lonsite <- lonsite - 360
      #print(paste0(site_list[s], " convert lon from degrees east to range [-180,180] | longitude is now: ", lonsite))
    }
    
    # create empty vectors based on the length of the 'swdown' meteo vector
    pdep=co2_dC13=co2_DC14=nhx=noy <- numeric(length(swdown))
    # loop over forcing years
    for (year in c(startyear:endyear)) {
      print(paste0(scenarioFileID, " ", year))
      
      ## N deposition
      #  this data covers the period 1850 to 2100
      if (year < 1850){ 
        fid<-nc_open("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_1850.nc")
      } else if (year < 1996){
        # use Ndep data created by Sönke
        fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_",year,".nc"))
      } else {
        fid<-nc_open(paste0("/Net/Groups/BGI/model/OCN/input/gridded/NDEP/CESM-CAM/ndep_RCP8.5_",year,".nc"))
      }
      lonndep   <- ncvar_get(fid,'lon')
      latndep   <- ncvar_get(fid,'lat')
      latx      <- which.min(abs(latndep-as.numeric(latsite)))
      lonx      <- which.min(abs(lonndep-as.numeric(lonsite)))
      nhx_dummy <- ncvar_get(fid,'NHx_deposition')[lonx,latx,]*86400*1000*1000
      noy_dummy <- ncvar_get(fid,'NOy_deposition')[lonx,latx,]*86400*1000*1000
      nhx_ratio <- nhx_dummy/(nhx_dummy+noy_dummy)
      nhx[1:(dpm[1]*48)+365*48*(year-startyear)] <- rep(nhx_ratio[1],dpm[1]*48)
      noy[1:(dpm[1]*48)+365*48*(year-startyear)] <- rep((1-nhx_ratio[1]),dpm[1]*48)
      for (m in 2:12) {
        nhx[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep(nhx_ratio[m],dpm[m]*48)
        noy[(1:(dpm[m]*48))+sum(dpm[1:(m-1)])*48+365*48*(year-startyear)]<-rep((1-nhx_ratio[m]),dpm[m]*48)
      }
      
      nc_close(fid)
      
      #CO2
      if(latsite>30) C14lat=2
      if(latsite>-30&latsite<=30) C14lat=3
      if(latsite<=-30) C14lat=4
      if (year < 1850){ 
        co2_dC13[1:17520+365*48*(year-startyear)] <- rep(co2y_dC13[(co2y_dC13[,1]-0.5)==1850,2],17520)
        co2_DC14[1:17520+365*48*(year-startyear)] <- rep(co2y_DC14[(co2y_DC14[,1]-0.5)==1850,C14lat],17520)
      } else if (year > 2019){ 
        co2_dC13[1:17520+365*48*(year-startyear)] <- rep(co2y_dC13[(co2y_dC13[,1]-0.5)==2019,2],17520)
        co2_DC14[1:17520+365*48*(year-startyear)] <- rep(co2y_DC14[(co2y_DC14[,1]-0.5)==2019,C14lat],17520)
      } else {
        co2_dC13[1:17520+365*48*(year-startyear)] <- rep(co2y_dC13[(co2y_dC13[,1]-0.5)==year,2],17520)
        co2_DC14[1:17520+365*48*(year-startyear)] <- rep(co2y_DC14[(co2y_DC14[,1]-0.5)==year,C14lat],17520)
      }
      
      ## P deposition 
      ## using a fixed value (usually from "/Net/Groups/BGI/work_2/QUINCY/SLM/InputDataSources/P-DEP/nitrogenandphosphorus2x2annualdep.nc")
      ## add extra P for years 2020:2022 for the two scenarios
      # MDP: Moderate P addition
      # HIP: High P addition
      # default is NOP
      phosphorus <- 9.3 / 365.0  # value is in 'P per day'
      # modify default P deposition for years 2020:2022 if MDP or HIP
      if (year %in% c(2020:2022)){
        if (grepl("MDP", scenario_list[s])){
          phosphorus <- 1500 / 365.0
          print(paste0("increase P deposition to 1500/365 (MDP scenario) for year: ", year))
        } else if (grepl("HIP", scenario_list[s])){
          phosphorus <- 3000 / 365.0
          print(paste0("increase P deposition to 3000/365 (HIP scenario) for year: ", year))
        }
      }
      pdep[1:17520+365*48*(year-startyear)] <- rep(phosphorus, 17520)
    }
    
    nhx <- nhx*Ndep*1000/365
    noy <- noy*Ndep*1000/365
    
    
    # combine individual variables into one array for output
    forcing<-cbind(year_out[1:nstep],doy_out[1:nstep],hod_out[1:nstep],swdown[1:nstep],lwdown[1:nstep],tair[1:nstep],qair[1:nstep],
                   press[1:nstep],rain[1:nstep],snow[1:nstep],wind[1:nstep],co2[1:nstep],co2_dC13[1:nstep],co2_DC14[1:nstep],
                   nhx[1:nstep],noy[1:nstep],pdep[1:nstep])
    
    # write output append yearly data (header created above)
    print(" writing file")
    write.table(format(round(forcing,digits=4),nsmall=1), file=outScenarioFile, append=TRUE, quote=FALSE,
                col.names=FALSE, row.names=FALSE, sep=" ")
    
  } else {
    print(paste0("site file not available: ", inSiteFile))
  }
}

print(" Done !")

## end ##
