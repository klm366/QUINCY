#!/usr/bin/env Rscript

## Rscript for creating figures comparing SLM simulated GPP & NPP of ALL CRUNCEP sites with global forest data base
##
## Uses input from:
##  1) observational data
##  2) the script   providing simulated data across sites
##    a) mean values:         slm_site_mean.csv
##    a) standard deviation:  slm_site_sd.csv
##
## author Sönke Zaehle @ MPI-BGC Jena
#############################################################################################################

## pass from bash the path to the simulation results main folder
args<-commandArgs(TRUE)
ipath<-args[1]  # /Net/Groups/BGI/scratch/""$username""/slm/"exp_ID"/
setwd(ipath)

## get observational data - absolute path
obs_path <- '/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/'
obs_file <- 'global_forest_database_all_sites_01082017.csv'
obs      <- read.csv(paste0(obs_path,obs_file), header=T, stringsAsFactors=F)
obs_nvar <- ncol(obs)

## get model results 
model_mean_file <- 'slm_site_mean.csv'
model_std_file  <- 'slm_site_sd.csv'
model_mean      <- read.csv(model_mean_file, header=T, stringsAsFactors=F)
model_std       <- read.csv(model_std_file, header=T, stringsAsFactors=F)
nsite_all       <- nrow(model_mean)
nsite_sel       <- 0                  # selected sites, present in 'obs' and 'model' data

## only sites present in 'obs' and 'model' data - potentially changing model_mean & model_std objects
print("  Note: this script does create figures only for GPP & NPP for ALL cruncep sites present in the the 'obs_file'.")

# get a vector with only the sites that are present in 'obs'
# and also exclude sites present in 'obs' but not in 'model_mean'
mod2obs_idx     <- numeric(nsite_all)
model_sel_sites <- numeric(nsite_all)
for (s in 1:nsite_all){
  if (length(which(model_mean$Site_ID[s] == obs$Site.Code)) > 0){
    nsite_sel                  <- nsite_sel + 1
    model_sel_sites[nsite_sel] <- s
    mod2obs_idx[s]             <- which(model_mean$Site_ID[s] == obs$Site.Code)
  }
}

# exit if no site was selected
if (nsite_sel == 0){
  print("  No simulated site matches the oberserved site list. EXIT this script.")
  quit(save="no", status=1)
}

model_sel_sites <- model_sel_sites[1:nsite_sel]
mod2obs_idx     <- mod2obs_idx[1:nsite_sel]
model_mean      <- model_mean[model_sel_sites,]
model_std       <- model_std[model_sel_sites,]

# create data.frame with obs data matching the model data
obs4mod     <- obs[1:nsite_sel,]   # create a new data.frame by copying part of the obs data.frame
for (s in mod2obs_idx){
  for (v in 1:obs_nvar){
    obs4mod[s,v]<-obs[mod2obs_idx[s],v]
  }
}

# write file for documentation
write.csv(obs4mod,file='obs_all_site.csv',row.names=F)


## prepare figures
pchidx<-c(1,19,0,15,2,17)
pchleg<-c('BroadDec','BroadEve','NeedleDec','NeedleEve','BroadMix','Mixed')
pch_site<-numeric(nsite_sel)
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Deciduous']<-pchidx[1]
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Evergreen']<-pchidx[2]
pch_site[obs4mod$Needles.leaves=='Needle-leaved'&obs4mod$Evergreen.Deciduous=='Deciduous']<-pchidx[3]
pch_site[obs4mod$Needles.leaves=='Needle-leaved'&obs4mod$Evergreen.Deciduous=='Evergreen']<-pchidx[4]
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Mixed']<-pchidx[5]
pch_site[obs4mod$Needles.leaves=='Mixed'&obs4mod$Evergreen.Deciduous=='Mixed']<-pchidx[6]

colidx<-c('blue','black','orange','red')
colleg<-c('boreal','temperate','mediterr.','tropical')
col_site<-numeric(nsite_sel)
col_site[obs4mod$Climatic.region=="Boreal Humid"|obs4mod$Climatic.region=="Boreal Semi-arid"]<-colidx[1]
col_site[obs4mod$Climatic.region=="Temperate Humid"|obs4mod$Climatic.region=="Temperate Semi-arid"]<-colidx[2]
col_site[obs4mod$Climatic.region=="Mediterranean Warm"]<-colidx[3]
col_site[obs4mod$Climatic.region=="Tropical Humid"|obs4mod$Climatic.region=="Tropical Semi-arid"]<-colidx[4]

plot_mod2obs<-function(xvar,yvar,label) {
  vmax<-max(xvar,yvar,na.rm=TRUE)
  vmax<-vmax+(vmax*0.2)  # add 20% to create space for the legend
  plot(xvar,yvar,xlim=c(0,vmax),ylim=c(0,vmax),pch=pch_site,col=col_site,
       xlab=paste('observed',label),ylab=paste('modelled',label))
  abline(0,1)
  legend('bottomright',legend=c(pchleg,colleg),pch=c(pchidx,rep(19,4)),col=c(rep('grey',6),colidx),cex=0.7)
}

plot_sitelocs<-function(lon,lat,pchx,colx) {

   library(maps)
   library(lattice)

   tatdegx=seq(-180,180,45)
   tatdegy=seq(-180,180,30)
   frgmap=map('world',xlim=c(-180,180),ylim=c(-60,90),plot=FALSE,interior=TRUE,add=TRUE)

   value=rep(0,10);lon1=value;lat1=value
   grid=data.frame(lon1,lat1,value);gridfmla=value~lon1*lat1

   pplot=levelplot(gridfmla,data=grid,region=TRUE,xlim=c(-180,180),
        ylim=c(-60,90),aspect=0.5,pretty=TRUE,
        colorkey=FALSE,col.regions='transparent',
        xlab='longitude',ylab='latitude',
        scales=list(x=list(at=tatdegx),y=list(at=tatdegy),tck=1,cex=1))

   print(pplot)

   panel.locs=trellis.currentLayout()
   for (i in 1:nrow(panel.locs)) {
         for (j in 1:ncol(panel.locs)) {
          if (panel.locs[i,j]>0) {
             trellis.focus("panel",j,i,highlight=FALSE);
             panel.lines(frgmap,col="grey",cex=0.1,lwd=1)
             for ( l in 1:length(lon) ) { 
                panel.points(lon[l],lat[l],pch=pchx[l],col=colx[l])
                 }
               }
         }
   }
   trellis.unfocus()
}


pdf('compare_sim_cruncep_all_vs_gfdb-plots.pdf')

plot_sitelocs(obs4mod$Longitude,obs4mod$Latitude,pch_site,col_site)

xvar<-obs4mod$TNPPcorr;yvar=model_mean$NPP*12.001;label='NPP (observ corrected) [gC/m2/yr]'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$TNPP;yvar=model_mean$NPP*12.001;label='NPP [gC/m2/yr]'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$GPP;yvar=model_mean$GPP*12.001;label='GPP [gC/m2/yr]'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$TNPPcorr;yvar=model_mean$NPP/model_mean$GPP;label='NPP/GPP (NPP observ corrected)' 
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$TNPP;yvar=model_mean$NPP/model_mean$GPP;label='NPP/GPP' 
plot_mod2obs(xvar,yvar,label)


trash <- dev.off()
