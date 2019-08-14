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

## Rscript for reading SLM results provided by simulations at the HPC via contrib/site_level/prep_run_slm.sh
##  calculating both mean and sd values and saving them in two individual files 
## These files will be used to compare SLM simulation results of CRUNCEP sites with global forest data base

## Option in prep_run_slm.sh: compareSimCRUNCEPvsGFDB="T"
## 
## author Sönke Zaehle @ MPI-BGC Jena
#############################################################################################################


## pass from bash the path to the simulation results main folder
args  <- commandArgs(TRUE)
ipath <- args[1]  # /Net/Groups/BGI/scratch/""$username""/slm/"exp_ID"/
setwd(ipath)


## detect sites
# get all directories within the results folder - relative paths
dirnames  <- system('ls -d */',intern=T)
# remove folders from list that contain other stuff than the simulation results
dirnames  <- dirnames[!grepl("x_binary",dirnames)]
dirnames  <- dirnames[!grepl("parameter_list",dirnames)]
dirnames  <- dirnames[!grepl("default_files",dirnames)]
dirnames  <- dirnames[!grepl("parameter_sensitivity_files",dirnames)]
dirnames  <- dirnames[!grepl("out",dirnames)]
dirnames  <- dirnames[!grepl("error",dirnames)]
# number of directories/sites
nsite     <- length(dirnames)
# create two arrays of length(directories and one column) that will be filled with the directory's names and Site-ID's (see loop)
site_str  <- array(NA,c(nsite,1))
colnames(site_str) <- 'Site_ID'

## loop over sites
for (s in 1:nsite) {
   print(paste(s,dirnames[[s]]))
   site_str[s]<-tail(unlist(strsplit(dirnames[[s]],'/')),n=1)

   cpool            <- read.table(paste(dirnames[[s]],'output/','vegpoolC_daily.txt',sep=''),header=T)
   c1314pool        <- read.table(paste(dirnames[[s]],'output/','vegpooldC13DC14_daily.txt',sep=''),header=T)
   npool            <- read.table(paste(dirnames[[s]],'output/','vegpoolN_daily.txt',sep=''),header=T)
   n15pool          <- read.table(paste(dirnames[[s]],'output/','vegpooldN15_daily.txt',sep=''),header=T)
   ppool            <- read.table(paste(dirnames[[s]],'output/','vegpoolP_daily.txt',sep=''),header=T)
   h2opool          <- read.table(paste(dirnames[[s]],'output/','soilpoolPhysics_daily.txt',sep=''),header=T)
   diags            <- read.table(paste(dirnames[[s]],'output/','veg_diagnostics_daily.txt',sep=''),header=T)
   cflux            <- read.table(paste(dirnames[[s]],'output/','vegfluxC_daily.txt',sep=''),header=T)
   nflux            <- read.table(paste(dirnames[[s]],'output/','vegfluxN_daily.txt',sep=''),header=T)
   pflux            <- read.table(paste(dirnames[[s]],'output/','vegfluxP_daily.txt',sep=''),header=T)
   cnp_soilpool     <- read.table(paste(dirnames[[s]],'output/','soilpoolOrg_daily.txt',sep=''),header=T)
   c13c14n15_soilpool <- read.table(paste(dirnames[[s]],'output/','soilpooldC13DC14dN15_daily.txt',sep=''),header=T)
   cc13c14_soilflux <- read.table(paste(dirnames[[s]],'output/','soilfluxCdC13DC14_daily.txt',sep=''),header=T)
   n_soilflux       <- read.table(paste(dirnames[[s]],'output/','soilfluxN_daily.txt',sep=''),header=T)  
   p_soilflux       <- read.table(paste(dirnames[[s]],'output/','soilfluxP_daily.txt',sep=''),header=T)  
   atm_forcing      <- read.table(paste(dirnames[[s]],'output/','atm_forcing_daily.txt',sep=''),header=T) 
   h2oflux          <- read.table(paste(dirnames[[s]],'output/','vegfluxH2O_daily.txt',sep=''),header=T)  

   ntime            <- dim(cpool)[1]
   ncpoolvar        <- dim(cpool)[2]-1
   ncfluxvar        <- dim(cflux)[2]-1
   nnfluxvar        <- dim(nflux)[2]-1
   npfluxvar        <- dim(pflux)[2]-1
   ncnpsoilpoolvar  <- dim(cnp_soilpool)[2]-1
   nc13c14n15soilpoolvar <- dim(c13c14n15_soilpool)[2]-1
   ncsoilfluxvar    <- dim(cc13c14_soilflux)[2]-1
   nnsoilfluxvar    <- dim(n_soilflux)[2]-1
   npsoilfluxvar    <- dim(p_soilflux)[2]-1
   natmforcingvar   <- dim(atm_forcing)[2]-1
   nh2ofluxvar      <- dim(h2oflux)[2]-1
   nyear            <- as.integer(round(ntime/365,0))
   yearidx          <- numeric(ntime)  # creates a vector of length ntime
   years            <- numeric(nyear)  # creates a vector of length nyear
   
   # add nyear sets of integer values of the particular year (usually 1 to 30)
   for (y in 1:nyear) {
     yearidx[(1:365)+(y-1)*365]<-y
   }
   # and correct for length of yearidx
   yearidx <- yearidx[1:ntime]

   ## default period of data used for analyses is 1:nyear (with nyear being usually 30)
   # you could change it to a different period within this default period, such as 10:19
   startYrDat <- 1
   lastYrDat  <- nyear
   ##

   # Cpools
   if(s==1) { 
      cpool_mean<-array(NA,c(nsite,ncpoolvar))
      colnames(cpool_mean)<-colnames(cpool)[1+1:ncpoolvar]
      cpool_std<-cpool_mean
   }
   for (var in 1:ncpoolvar) {
     hlp1<-aggregate(cpool[,var+1],list(yearidx),max) 
     cpool_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     cpool_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   # Cfluxes
   if(s==1) { 
      cflux_mean<-array(NA,c(nsite,ncfluxvar))
      colnames(cflux_mean)<-colnames(cflux)[1+1:ncfluxvar]
      cflux_std<-cflux_mean
   }
   for (var in 1:ncfluxvar) {
     hlp1<-aggregate(cflux[,var+1],list(yearidx),sum) 
     cflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     cflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   # Nfluxes
   if(s==1) { 
      nflux_mean<-array(NA,c(nsite,nnfluxvar))
      colnames(nflux_mean)<-colnames(nflux)[1+1:nnfluxvar]
      nflux_std<-nflux_mean
   }
   for (var in 1:nnfluxvar) {
     hlp1<-aggregate(nflux[,var+1],list(yearidx),sum) 
     nflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     nflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   # Pfluxes
   if(s==1) { 
      pflux_mean<-array(NA,c(nsite,npfluxvar))
      colnames(pflux_mean)<-colnames(pflux)[1+1:npfluxvar]
      pflux_std<-pflux_mean
   }
   for (var in 1:npfluxvar) {
     hlp1<-aggregate(pflux[,var+1],list(yearidx),sum) 
     pflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     pflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   
   #leaf & labile CNP
   if(s==1) { 
      stoch_mean<-array(NA,c(nsite,6))
      colnames(stoch_mean)<-c('Leaf_CN','Labile_CN','Reserve_CN','Leaf_NP','Labile_NP','Reserve_NP')
      stoch_std<-stoch_mean
   }
   var_idx<-which(colnames(cpool)=='Leaf_C'|colnames(cpool)=='Labile_C'|colnames(cpool)=='Reserve_C')
   res_idx<-1
   for (var in var_idx) {
     hlp1<-aggregate(cpool[,var],list(yearidx),sum) 
     hlp2<-aggregate(npool[,var],list(yearidx),sum) 
     hlp3<-aggregate(ppool[,var],list(yearidx),sum) 
     stoch_mean[s,res_idx]<-mean(hlp1$x[startYrDat:lastYrDat]/hlp2$x[startYrDat:lastYrDat],na.rm=TRUE)
     stoch_mean[s,res_idx+3]<-mean(hlp2$x[startYrDat:lastYrDat]/hlp3$x[startYrDat:lastYrDat],na.rm=TRUE)
     stoch_std[s,res_idx]<-sd(hlp1$x[startYrDat:lastYrDat]/hlp2$x[startYrDat:lastYrDat],na.rm=TRUE)
     stoch_std[s,res_idx+3]<-sd(hlp2$x[startYrDat:lastYrDat]/hlp3$x[startYrDat:lastYrDat],na.rm=TRUE)
     res_idx<-res_idx+1
   }

   #isotope signature
   if(s==1) { 
      c_isotope_mean<-array(NA,c(nsite,3))
      colnames(c_isotope_mean)<-c('Leaf_dC13','SapWood_DC14','HeartWood_DC14')
      c_isotope_std<-c_isotope_mean
   }
   var_idx<-which(colnames(c1314pool)=='Leaf_dC13'|colnames(c1314pool)=='SapWood_DC14'|colnames(c1314pool)=='HeartWood_DC14')
   for (var in var_idx) {
     cn<-colnames(c1314pool)
     if(cn[var]=='Leaf_dC13') c1314pool[abs(c1314pool[,var])<0.0001,var]<-NA
     hlp1<-aggregate(c1314pool[,var],list(yearidx),mean,na.rm=T)
     c_isotope_mean[s,which(colnames(c_isotope_mean)==colnames(c1314pool)[var])]<-mean(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
     c_isotope_std[s,which(colnames(c_isotope_mean)==colnames(c1314pool)[var])]<-sd(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
   }

   if(s==1) { 
      n_isotope_mean<-array(NA,c(nsite,2))
      colnames(n_isotope_mean)<-c('Leaf_dN15','Total_dN15')
      n_isotope_std<-n_isotope_mean
   }
   var_idx<-which(colnames(n15pool)=='Leaf_dN15'|colnames(n15pool)=='Total_dN15')
   for (var in var_idx) {
     cn<-colnames(n15pool)
     if(cn[var]=='Leaf_dN15') n15pool[abs(n15pool[,var])<0.0001,var]<-NA
     hlp1<-aggregate(n15pool[,var],list(yearidx),mean,na.rm=T)
     n_isotope_mean[s,which(colnames(n_isotope_mean)==colnames(n15pool)[var])]<-mean(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
     n_isotope_std[s,which(colnames(n_isotope_mean)==colnames(n15pool)[var])]<-sd(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
   }

   #soil water
   if(s==1) { 
      soilwater_mean<-array(NA,c(nsite,2))
      colnames(soilwater_mean)<-c('RootZoneSoilMoisture','RootZoneSoilWaterPotential')
      soilwater_std<-soilwater_mean
   }
   var_idx<-which(colnames(h2opool)=='RootZoneSoilMoisture'|colnames(h2opool)=='RootZoneSoilWaterPotential')
   for (var in var_idx) {
     if(colnames(h2opool)[var]=='RootZoneSoilMoisture') {
         maxval<-max(h2opool[,var])
         tsval<-h2opool[,var]/maxval
         hlp1<-aggregate(tsval,list(yearidx),mean)
     }
     else {hlp1<-aggregate(h2opool[,var],list(yearidx),min)} 
     soilwater_mean[s,which(colnames(soilwater_mean)==colnames(h2opool)[var])]<-mean(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
     soilwater_std[s,which(colnames(soilwater_mean)==colnames(h2opool)[var])]<-sd(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
   }

   #other diagnostics
   if(s==1) { 
      diag_mean<-array(NA,c(nsite,3))
      colnames(diag_mean)<-c('LAI','Height','ind_density')
      diag_std<-diag_mean
   }
   var_idx<-which(colnames(diags)=='LAI'|colnames(diags)=='Height'|colnames(diags)=='ind_density')
   for (var in var_idx) {
     if(colnames(diags)[var]=='LAI') {hlp1<-aggregate(diags[,var],list(yearidx),max)} 
     else {hlp1<-aggregate(diags[,var],list(yearidx),mean)} 
     diag_mean[s,which(colnames(diag_mean)==colnames(diags)[var])]<-mean(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
     diag_std[s,which(colnames(diag_mean)==colnames(diags)[var])]<-sd(hlp1$x[startYrDat:lastYrDat],na.rm=TRUE)
   }

   # soil organic pools
   if(s==1) { 
      cnp_soilpool_mean<-array(NA,c(nsite,ncnpsoilpoolvar))
      colnames(cnp_soilpool_mean)<-colnames(cnp_soilpool)[1+1:ncnpsoilpoolvar]
      cnp_soilpool_std<-cnp_soilpool_mean
   }
   for (var in 1:ncnpsoilpoolvar) {
     hlp1<-aggregate(cnp_soilpool[,var+1],list(yearidx),max) 
     cnp_soilpool_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     cnp_soilpool_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   # isotopic signature soil organic pools
   if(s==1) { 
      c13c14n15_soilpool_mean<-array(NA,c(nsite,nc13c14n15soilpoolvar))
      colnames(c13c14n15_soilpool_mean)<-colnames(c13c14n15_soilpool)[1+1:nc13c14n15soilpoolvar]
      c13c14n15_soilpool_std<-c13c14n15_soilpool_mean
   }
   for (var in 1:nc13c14n15soilpoolvar) {
     hlp1<-aggregate(c13c14n15_soilpool[,var+1],list(yearidx),max) 
     c13c14n15_soilpool_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     c13c14n15_soilpool_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
   # soil cfluxes pools
   if(s==1) { 
      cc13c14_soilflux_mean<-array(NA,c(nsite,ncsoilfluxvar))
      colnames(cc13c14_soilflux_mean)<-colnames(cc13c14_soilflux)[1+1:ncsoilfluxvar]
      cc13c14_soilflux_std<-cc13c14_soilflux_mean
   }
   for (var in 1:ncsoilfluxvar) {
     hlp1<-aggregate(cc13c14_soilflux[,var+1],list(yearidx),sum) 
     cc13c14_soilflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     cc13c14_soilflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }

   # soil n fluxes
   if(s==1) { 
      n_soilflux_mean<-array(NA,c(nsite,nnsoilfluxvar))
      colnames(n_soilflux_mean)<-colnames(n_soilflux)[1+1:nnsoilfluxvar]
      n_soilflux_std<-n_soilflux_mean
   }
   for (var in 1:nnsoilfluxvar) {
     hlp1<-aggregate(n_soilflux[,var+1],list(yearidx),sum) 
     n_soilflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     n_soilflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }

   # soil p fluxes
   if(s==1) { 
      p_soilflux_mean<-array(NA,c(nsite,npsoilfluxvar))
      colnames(p_soilflux_mean)<-colnames(p_soilflux)[1+1:npsoilfluxvar]
      p_soilflux_std<-p_soilflux_mean
   }
   for (var in 1:npsoilfluxvar) {
     hlp1<-aggregate(p_soilflux[,var+1],list(yearidx),sum) 
     p_soilflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     p_soilflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }

   # atmospheric forcing 
   if(s==1) { 
      atm_forcing_mean<-array(NA,c(nsite,natmforcingvar))
      colnames(atm_forcing_mean)<-colnames(atm_forcing)[1+1:natmforcingvar]
      atm_forcing_std<-atm_forcing_mean
   }
   for (var in 1:natmforcingvar) {
     hlp1<-aggregate(atm_forcing[,var+1],list(yearidx),mean) 
     atm_forcing_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     atm_forcing_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }

   # h2o fluxes 
   if(s==1) { 
      h2oflux_mean<-array(NA,c(nsite,nh2ofluxvar))
      colnames(h2oflux_mean)<-colnames(h2oflux)[1+1:nh2ofluxvar]
      h2oflux_std<-h2oflux_mean
   }
   for (var in 1:nh2ofluxvar) {
     hlp1<-aggregate(h2oflux[,var+1],list(yearidx),sum) 
     h2oflux_mean[s,var]<-mean(hlp1$x[startYrDat:lastYrDat])
     h2oflux_std[s,var]<-sd(hlp1$x[startYrDat:lastYrDat])
   }
}   

# hack for simple
soil_org_isotope_mean<-array(NA,c(nsite,6))
colnames(soil_org_isotope_mean)<-c('Litter_dC13','Litter_DC14','Litter_dN15','SOM_dC13','SOM_DC14','SOM_dN15')
soil_org_isotope_std<-soil_org_isotope_mean

pool_mean<-data.frame(cnp_soilpool_mean)
pool_std<-data.frame(cnp_soilpool_std)
iso_mean<-data.frame(c13c14n15_soilpool_mean)
iso_std<-data.frame(c13c14n15_soilpool_std)
soil_org_isotope_mean[,1]<-((iso_mean$MetLitterdC13*pool_mean$MetLitterC+
                             iso_mean$StrLitterdC13*pool_mean$StrLitterC+
                             iso_mean$WodLitterdC13*pool_mean$WodLitterC)/
                            (pool_mean$MetLitterC+pool_mean$StrLitterC+pool_mean$WodLitterC))
soil_org_isotope_mean[,2]<-((iso_mean$MetLitterDC14*pool_mean$MetLitterC+
                             iso_mean$StrLitterDC14*pool_mean$StrLitterC+
                             iso_mean$WodLitterDC14*pool_mean$WodLitterC)/
                            (pool_mean$MetLitterC+pool_mean$StrLitterC+pool_mean$WodLitterC))
soil_org_isotope_mean[,3]<-((iso_mean$MetLitterdN15*pool_mean$MetLitterN+
                             iso_mean$StrLitterdN15*pool_mean$StrLitterN+
                             iso_mean$WodLitterN15*pool_mean$WodLitterN)/
                            (pool_mean$MetLitterN+pool_mean$StrLitterN+pool_mean$WodLitterN))
soil_org_isotope_mean[,4]<-((iso_mean$FastdC13*pool_mean$FastC+
                             iso_mean$SlowdC13*pool_mean$SlowC)/
                            (pool_mean$FastC+pool_mean$SlowC))
soil_org_isotope_mean[,5]<-((iso_mean$FastDC14*pool_mean$FastC+
                             iso_mean$SlowDC14*pool_mean$SlowC)/
                            (pool_mean$FastC+pool_mean$SlowC))
soil_org_isotope_mean[,6]<-((iso_mean$FastdN15*pool_mean$FastN+
                             iso_mean$SlowdN15*pool_mean$SlowN)/
                            (pool_mean$FastN+pool_mean$SlowN))


# mean values:         slm_site_mean.csv
out_mean<-cbind(site_str,cpool_mean,stoch_mean,diag_mean,cnp_soilpool_mean,
                         cflux_mean,nflux_mean,n_soilflux_mean,pflux_mean,p_soilflux_mean,h2oflux_mean,
                         c_isotope_mean,n_isotope_mean,soil_org_isotope_mean,cc13c14_soilflux_mean,
                         soilwater_mean,atm_forcing_mean)
write.csv(out_mean,file='slm_site_mean.csv',quote=F)
# standard deviation:  slm_site_sd.csv
out_std<-cbind(site_str,cpool_std,stoch_std,diag_std,cnp_soilpool_std,
                         cflux_std,nflux_std,n_soilflux_std,pflux_std,p_soilflux_std,h2oflux_std,
                         c_isotope_std,n_isotope_std,soil_org_isotope_std,cc13c14_soilflux_std,
                         soilwater_std,atm_forcing_std)
write.csv(out_std,file='slm_site_sd.csv',quote=F)
