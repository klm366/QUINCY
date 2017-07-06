#!/usr/bin/env Rscript

## Rscript for reading SLM results provided by simulations at the HPC via contrib/site_level/prep_run_slm.sh
##  calculating both mean and sd values and saving them in two individual files 
## These files will be used to compare SLM simulation results of CRUNCEP sites with global forest data base

## Option in prep_run_slm.sh: compareSimCRUNCEPvsGFDB="T"
## 
## author Sönke Zaehle @ MPI-BGC Jena
#############################################################################################################


## pass from bash the path to the simulation results main folder
args<-commandArgs(TRUE)
ipath<-args[1]  # /Net/Groups/BGI/scratch/""$username""/slm/"exp_ID"/
setwd(ipath)


## detect sites
# get all directories within the results folder - relative paths
dirnames<-system('ls -d */',intern=T)
# remove folders from list that contain other stuff than the simulation results
dirnames<-dirnames[!grepl("x_binary",dirnames)]
dirnames<-dirnames[!grepl("parameter_list",dirnames)]
dirnames<-dirnames[!grepl("default_files",dirnames)]
dirnames<-dirnames[!grepl("parameter_sensitivity_files",dirnames)]
dirnames<-dirnames[!grepl("out",dirnames)]
dirnames<-dirnames[!grepl("error",dirnames)]
# number of directories/sites
nsite<-length(dirnames)
# create two arrays of length(directories and one column) that will be filled with the directory's names and Site-ID's (see loop)
site_str<-array(NA,c(nsite,1))
colnames(site_str)<-'Site_ID'

## loop over sites
for (s in 1:nsite) {
   print(paste(s,dirnames[[s]]))
   site_str[s]<-tail(unlist(strsplit(dirnames[[s]],'/')),n=1)

   cpool<-read.table(paste(dirnames[[s]],'output/','vegpoolC_daily.txt',sep=''),header=T)
   c1314pool<-read.table(paste(dirnames[[s]],'output/','vegpooldC13dC14_daily.txt',sep=''),header=T)
   npool<-read.table(paste(dirnames[[s]],'output/','vegpoolN_daily.txt',sep=''),header=T)
   ppool<-read.table(paste(dirnames[[s]],'output/','vegpoolP_daily.txt',sep=''),header=T)
   h2opool<-read.table(paste(dirnames[[s]],'output/','vegpoolH2O_daily.txt',sep=''),header=T)
   diags<-read.table(paste(dirnames[[s]],'output/','veg_diagnostics_daily.txt',sep=''),header=T)
   cflux<-read.table(paste(dirnames[[s]],'output/','vegfluxC_daily.txt',sep=''),header=T)
   nflux<-read.table(paste(dirnames[[s]],'output/','vegfluxN_daily.txt',sep=''),header=T)
   pflux<-read.table(paste(dirnames[[s]],'output/','vegfluxP_daily.txt',sep=''),header=T)
   cnp_soilpool<-read.table(paste(dirnames[[s]],'output/','soilpoolOrg_daily.txt',sep=''),header=T)
   c13c14_soilpool<-read.table(paste(dirnames[[s]],'output/','soilpooldCdN_daily.txt',sep=''),header=T)
   cc13c14_soilflux<-read.table(paste(dirnames[[s]],'output/','soilfluxCdC13DC14_daily.txt',sep=''),header=T)
   np_soilflux<-read.table(paste(dirnames[[s]],'output/','soilfluxNP_daily.txt',sep=''),header=T)  
 
   ntime<-dim(cpool)[1]
   ncpoolvar<-dim(cpool)[2]-1
   ncfluxvar<-dim(cflux)[2]-1
   nnfluxvar<-dim(nflux)[2]-1
   npfluxvar<-dim(pflux)[2]-1
   ncnpsoilpoolvar<-dim(cnp_soilpool)[2]-1
   ncsoilfluxvar<-dim(cc13c14_soilflux)[2]-1
   nnpsoilfluxvar<-dim(np_soilflux)[2]-1
   nyear=ntime/365
   yearidx<- numeric(ntime)
   years<-numeric(nyear)
   for (y in 1:nyear) {
     yearidx[(1:365)+(y-1)*365]<-y
   }
   
   # Cpools
   if(s==1) { 
      cpool_mean<-array(NA,c(nsite,ncpoolvar))
      colnames(cpool_mean)<-colnames(cpool)[1+1:ncpoolvar]
      cpool_std<-cpool_mean
   }
   for (var in 1:ncpoolvar) {
     hlp1<-aggregate(cpool[,var+1],list(yearidx),max) 
     cpool_mean[s,var]<-mean(hlp1$x)
     cpool_std[s,var]<-sd(hlp1$x)
   }
   # Cfluxes
   if(s==1) { 
      cflux_mean<-array(NA,c(nsite,ncfluxvar))
      colnames(cflux_mean)<-colnames(cflux)[1+1:ncfluxvar]
      cflux_std<-cflux_mean
   }
   for (var in 1:ncfluxvar) {
     hlp1<-aggregate(cflux[,var+1],list(yearidx),sum) 
     cflux_mean[s,var]<-mean(hlp1$x)
     cflux_std[s,var]<-sd(hlp1$x)
   }
   # Nfluxes
   if(s==1) { 
      nflux_mean<-array(NA,c(nsite,nnfluxvar))
      colnames(nflux_mean)<-colnames(nflux)[1+1:nnfluxvar]
      nflux_std<-nflux_mean
   }
   for (var in 1:nnfluxvar) {
     hlp1<-aggregate(nflux[,var+1],list(yearidx),sum) 
     nflux_mean[s,var]<-mean(hlp1$x)
     nflux_std[s,var]<-sd(hlp1$x)
   }
   # Pfluxes
   if(s==1) { 
      pflux_mean<-array(NA,c(nsite,npfluxvar))
      colnames(pflux_mean)<-colnames(pflux)[1+1:npfluxvar]
      pflux_std<-pflux_mean
   }
   for (var in 1:npfluxvar) {
     hlp1<-aggregate(pflux[,var+1],list(yearidx),sum) 
     pflux_mean[s,var]<-mean(hlp1$x)
     pflux_std[s,var]<-sd(hlp1$x)
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
     stoch_mean[s,res_idx]<-mean(hlp1$x/hlp2$x,na.rm=TRUE)
     stoch_mean[s,res_idx+3]<-mean(hlp2$x/hlp3$x,na.rm=TRUE)
     stoch_std[s,res_idx]<-sd(hlp1$x/hlp2$x,na.rm=TRUE)
     stoch_std[s,res_idx+3]<-sd(hlp2$x/hlp3$x,na.rm=TRUE)
     res_idx<-res_idx+1
   }

   #isotope signature
   if(s==1) { 
      isotope_mean<-array(NA,c(nsite,3))
      colnames(isotope_mean)<-c('Labile_dC13','Labile_DC14','HeartWood_DC14')
      isotope_std<-isotope_mean
   }
   var_idx<-which(colnames(c1314pool)=='Labile_dC13'|colnames(c1314pool)=='Labile_DC14'|colnames(c1314pool)=='HeartWood_DC14')
   for (var in var_idx) {
     hlp1<-aggregate(c1314pool[,var],list(yearidx),mean)
     isotope_mean[s,which(colnames(isotope_mean)==colnames(c1314pool)[var])]<-mean(hlp1$x,na.rm=TRUE)
     isotope_std[s,which(colnames(isotope_mean)==colnames(c1314pool)[var])]<-sd(hlp1$x,na.rm=TRUE)
   }

   #soil water
   if(s==1) { 
      soilwater_mean<-array(NA,c(nsite,2))
      colnames(soilwater_mean)<-c('SoilMoisture','SoilWaterPotential')
      soilwater_std<-soilwater_mean
   }
   var_idx<-which(colnames(h2opool)=='SoilMoisture'|colnames(h2opool)=='SoilWaterPotential')
   for (var in var_idx) {
     if(colnames(h2opool)[var]=='SoilMoisture') {hlp1<-aggregate(h2opool[,var],list(yearidx),mean)/aggregate(h2opool[,var],list(yearidx),max)} 
     else {hlp1<-aggregate(h2opool[,var],list(yearidx),min)} 
     soilwater_mean[s,which(colnames(soilwater_mean)==colnames(h2opool)[var])]<-mean(hlp1$x,na.rm=TRUE)
     soilwater_std[s,which(colnames(soilwater_mean)==colnames(h2opool)[var])]<-sd(hlp1$x,na.rm=TRUE)
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
     diag_mean[s,which(colnames(diag_mean)==colnames(diags)[var])]<-mean(hlp1$x,na.rm=TRUE)
     diag_std[s,which(colnames(diag_mean)==colnames(diags)[var])]<-sd(hlp1$x,na.rm=TRUE)
   }

   # soil organic pools
   if(s==1) { 
      cnp_soilpool_mean<-array(NA,c(nsite,ncnpsoilpoolvar))
      colnames(cnp_soilpool_mean)<-colnames(cnp_soilpool)[1+1:ncnpsoilpoolvar]
      cnp_soilpool_std<-cnp_soilpool_mean
   }
   for (var in 1:ncnpsoilpoolvar) {
     hlp1<-aggregate(cnp_soilpool[,var+1],list(yearidx),max) 
     cnp_soilpool_mean[s,var]<-mean(hlp1$x)
     cnp_soilpool_std[s,var]<-sd(hlp1$x)
   }
   # soil cfluxes pools
   if(s==1) { 
      cc13c14_soilflux_mean<-array(NA,c(nsite,ncsoilfluxvar))
      colnames(cc13c14_soilflux_mean)<-colnames(cc13c14_soilflux)[1+1:ncsoilfluxvar]
      cc13c14_soilflux_std<-cc13c14_soilflux_mean
   }
   for (var in 1:ncsoilfluxvar) {
     hlp1<-aggregate(cc13c14_soilflux[,var+1],list(yearidx),sum) 
     cc13c14_soilflux_mean[s,var]<-mean(hlp1$x)
     cc13c14_soilflux_std[s,var]<-sd(hlp1$x)
   }
   # soil cfluxes pools
   if(s==1) { 
      np_soilflux_mean<-array(NA,c(nsite,nnpsoilfluxvar))
      colnames(np_soilflux_mean)<-colnames(np_soilflux)[1+1:nnpsoilfluxvar]
      np_soilflux_std<-np_soilflux_mean
   }
   for (var in 1:nnpsoilfluxvar) {
     hlp1<-aggregate(np_soilflux[,var+1],list(yearidx),sum) 
     np_soilflux_mean[s,var]<-mean(hlp1$x)
     np_soilflux_std[s,var]<-sd(hlp1$x)
   }
}   

# mean values:         slm_site_mean.csv
# standard deviation:  slm_site_sd.csv
out_mean<-cbind(site_str,cpool_mean,diag_mean,stoch_mean,cflux_mean,nflux_mean,pflux_mean,isotope_mean,soilwater_mean,cnp_soilpool_mean,cc13c14_soilflux_mean,np_soilflux_mean)
write.csv(out_mean,file='slm_site_mean.csv',quote=F)
out_std<-cbind(site_str,cpool_std,diag_std,stoch_std,cflux_std,nflux_std,pflux_std,isotope_std,soilwater_std,cnp_soilpool_std,cc13c14_soilflux_std,np_soilflux_std)
write.csv(out_std,file='slm_site_sd.csv',quote=F)
