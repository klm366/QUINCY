#!/usr/bin/env Rscript
# Interface 
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 3) {
  message("slm2ncdf.R requires longitude latitude and startyear as input. Terminating");quit()
}
lon <- as.numeric(Args[1])
lat <- as.numeric(Args[2])
year <- as.numeric(Args[3])

library(ncdf)

#function to creat cdf compliant files
create_cdf<-function(year,lon,lat,timeout,timeindex,fdata,miss_val,nvar,
                      fname,flongname,funit,fprec,fname_cdf) {

  #create dimensions
  dimLon<-dim.def.ncdf("longitude","degrees_east",lon,unlim=FALSE,create_dimvar=TRUE)
  dimLat<-dim.def.ncdf("latitude","degrees_north",lat,unlim=FALSE,create_dimvar=TRUE)
  dimTime<-dim.def.ncdf("time",paste("seconds since ",year[timeindex][1],"-01-01 00:00:00",sep=""),
                      timeout,unlim=TRUE,create_dimvar=TRUE)

  #create variables
  forvarid<-list()
  for (v in 1:nvar) {
    forvarid[[v]]<-var.def.ncdf(fname[v],funit[v],list(dimLon,dimLat,dimTime),
                                missval=miss_val,longname=flongname[v],prec=fprec[v])
  }

  #create and write file
  ncid<-create.ncdf(fname_cdf,forvarid)
  for (v in 1:nvar) {
     put.var.ncdf(ncid,forvarid[[v]],fdata[v,timeindex])
     att.put.ncdf(ncid,forvarid[[v]],"associate","time latitude longitude",prec="text")
     att.put.ncdf(ncid,forvarid[[v]],"axis","TYX",prec="text")
  }
  att.put.ncdf(ncid,0,"description","SLM derived FLUXNET simulation",prec='text')
  att.put.ncdf(ncid,0,"conventions","CDF",prec='text')
  att.put.ncdf(ncid,0,"date",date(),prec='text')
  close.ncdf(ncid)

}


fluxCfile<-read.table('vegfluxC_tile1_timestep.txt',header=T)
fluxH2Ofile<-read.table('vegfluxH2O_tile1_timestep.txt',header=T)
fluxEfile<-read.table('vegfluxEnergy_tile1_timestep.txt',header=T)

miss_val<--9.33e30
ntime<-length(fluxCfile[,1])
dt=round((fluxCfile[3,1]-fluxCfile[2,1])*86400)

variable_name <- array(c("PPFD","Rn","Qh","Qle","Ga","Gc","GPP","Reco","NEE","WUE","FAPAR"));nvar<-length(variable_name)
outname=outlongname=outunit=outprec=numeric(nvar)
outprec[]<-"double"
outvar<-array(0,c(nvar,ntime))
tstepaxis<-seq(1,ntime,1)
timeaxis<-tstepaxis*dt

fname<-'slm4fluxnet.nc'

for (v in 1:nvar) {
  outname[v]<-variable_name[v]
  if(variable_name[v]=='PPFD') {
     outlongname[v]<-"Photosynthetic Photon Flux Density"
     outunit[v]<-"µmol / m2 /s"
     outvar[v,]<-fluxEfile$PPFD
  }
  if(variable_name[v]=='Rn') {
     outlongname[v]<-"Net Radiation"
     outunit[v]<-"W / m2"
     outvar[v,]<-fluxEfile$Rnet
  }
  if(variable_name[v]=='Qh') {
     outlongname[v]<-"Sensible Heat FLux"
     outunit[v]<-"W / m2"
     outvar[v,]<-fluxEfile$Qh
  }
  if(variable_name[v]=='Qle') {
     outlongname[v]<-"Latent Heat FLux"
     outunit[v]<-"W / m2"
     outvar[v,]<-fluxEfile$Qle
  }
  if(variable_name[v]=='Ga') {
     outlongname[v]<-"Aerodynamic Conductance"
     outunit[v]<-"m/s"
     outvar[v,]<-fluxH2Ofile$GA
  }
  if(variable_name[v]=='Gc') {
     outlongname[v]<-"Canopy Conductance"
     outunit[v]<-"m/s"
     outvar[v,]<-fluxH2Ofile$GC
  }
  if(variable_name[v]=='GPP') {
     outlongname[v]<-"Gross Primary Production"
     outunit[v]<-"µmol C / m2 / s"
     outvar[v,]<-fluxCfile$GPP
  }
  if(variable_name[v]=='Reco') {
     outlongname[v]<-"Ecosystem Respiration"
     outunit[v]<-"µmol C / m2 / s"
     outvar[v,]<-fluxCfile$GrowthResp+fluxCfile$MaintResp
  }
  if(variable_name[v]=='NEE') {
     outlongname[v]<-"Net Ecosystem Production"
     outunit[v]<-"µmol C / m2 / s"
     outvar[v,]<--1*fluxCfile$NPP
  }
  if(variable_name[v]=='WUE') {
     outlongname[v]<-"Water-Use Efficiency"
     outunit[v]<-"gC / kg H20"
     outvar[v,]<-fluxCfile$GPP/fluxH2Ofile$Transpiration*12/1000000
     outvar[v,fluxH2Ofile$Transpiration<0.00000001]<-0
     outvar[v,is.na(outvar[v,])]<-0
  }
  if(variable_name[v]=='FAPAR') {
     outlongname[v]<-"Fraction of Absorbed Photosynthetically Active Radiation"
     outunit[v]<-"-"
     outvar[v,]<-fluxEfile$APAR/fluxEfile$PPFD
     outvar[v,is.na(outvar[v,])]<-0
  }

}

create_cdf(year,lon,lat,timeaxis,,outvar,miss_val,nvar,outname,outlongname,outunit,outprec,fname)


