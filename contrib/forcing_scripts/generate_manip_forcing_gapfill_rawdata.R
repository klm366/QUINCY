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
## preprocess the Manip / Majadas raw 1/2 h met forcing data provided by Tarek
## this raw data is gap filled, and it is ensured to have the correct number of data points (365 days per year)
##
## the original version of this script was written by Sönke (BGI/people/szaehle/Projects/QUINCY/GenerateMANP_Forcing_2018)
## 
## run this script by:
##  cd contrib/forcing_scripts/
##  Rscript --vanilla ./generate_manip_forcing_gapfill_rawdata.R
## 
## the output file of his script is needed by the generate_manip_static_slm_forcing.R
##
## output may be generated per default in the input folder where the raw data is located (aka workDir)
## 
#########################################################################################


## input (Rdata holds a R data.frame plus some meta-data)
#inRawData <- "ES-LMa_2009_2017.Rdata"
inRawData <- "ES-LMa_2009_2019.Rdata"

## output
#outGapFilledData <- "ES_LMa_2009_2017_gapfilled.Rdata"
outGapFilledData <- "ES_LMa_2009_2019_gapfilled.Rdata"

## working directory
workDir <- "/Net/Groups/BGI/work_2/QUINCY/SLM/Data/manip_majadas/raw_and_gapfilled_data/"

## message
print(paste0("Running gapfill on Manip/Majadas raw data: ",inRawData))

## load the dataFrame and store as df
load(paste0(workDir, inRawData))

## number of input years
nYears <- length( unique(df[,1]) )


## gapfill data - not for the last 5 columns which were just added by Tarek for QUINCY convenience, but without data
## (this data will be completed by the generate_manip_static_slm_forcing.R)
nvar<-length(df[1,])-5
dfnew<-df

for (v in 4:nvar) {
  if(any(is.na(df[,v]))) {
    miss<-which(is.na(df[,v]))
    for (m in 1:length(miss)) {
      dfnew[miss[m],v]<-mean(df[df$doy==df$doy[miss[m]]&df$hour==df$hour[miss[m]],v],na.rm=T)
    } 
  }
}


## remove 31 Dez for leap years (QUINCY does not consider leap years, for now [20191016])
# remove data for any 366th day
dfnew2<-dfnew[dfnew[,2]<366,]

## add data for the timesteps potentially missing in the last year of the raw data (needed for e.g. 2009-2017 dataset)
## create new data frame from dfnew2 with the number of rows as needed for all years to be complete 
## (assuming the last year of the raw data contains incomplete data not covering all days)
dfnew3<-as.data.frame(array(NA,c(nYears*365*48,17)))   # nYears * days * 48 (1/2 h timestep)
colnames(dfnew3)<-colnames(dfnew2)
dfnew3[1:length(dfnew2[,1]),]<-dfnew2
# calc the number of timesteps (1/2 h) missing
nmiss<-nrow(dfnew3) - nrow(dfnew2)

## if nmiss > 0
if (nmiss > 0){
  # get the sequence of days of the year that are missing
  doy_miss<-seq(365-nmiss/48+1,365,1)
  
  ## add "days" and "hours" information for the missing data
  d<-1
  h<-0.0
  for (i in 1:nmiss) {
     dfnew3[nrow(dfnew2) +i, 2] <- doy_miss[d]
     dfnew3[nrow(dfnew2) +i, 3] <- h
     h<-h+0.5
     if(h==24) {
       d<-d+1
       h<-0
     }
  }
  
  ## add missing data for the last year
  for (v in 4:nvar) {
  for (i in 1:nmiss) {
      dfnew3[length(dfnew2[,1])+i,v]<-mean(dfnew3[dfnew3$doy==dfnew3$doy[length(dfnew2[,1])+i]&dfnew3$hour==dfnew3$hour[length(dfnew2[,1])+i],v],na.rm=T)
    }
  }
  lastYear <- df[nrow(df),1]
  dfnew3[(length(dfnew2[,1])+1):length(dfnew3[,1]),1] <- lastYear
  

  # move the gapfilled data to 'df' object, as the load() in generate_manip_static_slm_forcing.R expects this object name
  df <- dfnew3
} else {  # if nmiss == 0
  # move the gapfilled data to 'df' object, as the load() in generate_manip_static_slm_forcing.R expects this object name
  df <- dfnew2
}

## correct the offset in the Majadas measures CO2 data in the initial about 2.5 years
# offset co2 is present from 1:42200 in Majadas data
print("CO2 offset correction for the initial ~ 2.5 years of measured data.")
aMean <- mean(df$co2_mixing_ratio[1:42200], na.rm = T)     # mean of the init ~ 2.5 years
bMean <- mean(df$co2_mixing_ratio[42201:59720], na.rm = T) # mean of the year after that
abOffset <- aMean - bMean  # it is a positive offset

## mean of global co2 for the years 2009-2012 from /Net/Groups/BGI/model/OCN/input/gridded/CO2/GCP_co2_global.dat
# 2009  385.80
# 2010  387.98
# 2011  389.79
# 2012  391.88
co2Global <- c(385.80, 387.98, 389.79, 391.88)
# global mean for the Majadas period [42201:59720] is assumed to be the mean of 2011/2012 global
globBMean <- mean (co2Global[3:4])
# correct the 2009-2011 (2011 only until 42200) and add the trend between the years to globBMean
rowSel <- which(df$year == 2009)
df$co2_mixing_ratio[rowSel] <- df$co2_mixing_ratio[rowSel] - abOffset - (globBMean - co2Global[1])
rowSel <- which(df$year == 2010)
df$co2_mixing_ratio[rowSel] <- df$co2_mixing_ratio[rowSel] - abOffset - (globBMean - co2Global[2])
rowSel <- which(df$year == 2011 )
rowSel <- rowSel[1:(42200-17520-17520)]
df$co2_mixing_ratio[rowSel] <- df$co2_mixing_ratio[rowSel] - abOffset - (globBMean - co2Global[3])


## write output
# this output file is used by the generate_manip_static_slm_forcing.R script to create actual QUINCY met-forcing
save(df, file=paste0(workDir,outGapFilledData))

## message
print("Done !")

# fin #
