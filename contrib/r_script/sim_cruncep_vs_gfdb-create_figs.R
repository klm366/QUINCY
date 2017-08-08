#!/usr/bin/env Rscript

## Rscript for creating figures comparing SLM simulation results of CRUNCEP sites with global forest data base
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
obs_path<-'/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/'
obs_file<-'global_forest_database_cue_sites_01082017.csv'
obs<-read.csv(paste0(obs_path,obs_file),stringsAsFactors=F)
obsdim<-dim(obs)

model_path<-''
model_expid<-''
model_mean_file<-'slm_site_mean.csv'
model_std_file<-'slm_site_sd.csv'
model_mean<-read.csv(paste(model_path,model_expid,model_mean_file,sep=''),stringsAsFactors=F)
model_std<-read.csv(paste(model_path,model_expid,model_std_file,sep=''),stringsAsFactors=F)

moddim<-dim(model_mean)
nsite<-moddim[1]

mod2obs_idx<-numeric(nsite)
obs4mod<-obs[1:nsite,]   # create a new data.frame by copying part of the obs data.frame
for (s in 1:nsite) {
  mod2obs_idx[s]<-which(model_mean$Site_ID[s]==obs$Site.Code)
  for (v in 1:obsdim[2]) obs4mod[s,v]<-obs[mod2obs_idx[s],v]
}
write.csv(obs4mod,file='obs_site.csv',row.names=F)

pchidx<-c(1,19,0,15,2,17)
pchleg<-c('BroadDec','BroadEve','NeedleDec','NeedleEve','BroadMix','Mixed')
pch_site<-numeric(nsite)
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Deciduous']<-pchidx[1]
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Evergreen']<-pchidx[2]
pch_site[obs4mod$Needles.leaves=='Needle-leaved'&obs4mod$Evergreen.Deciduous=='Deciduous']<-pchidx[3]
pch_site[obs4mod$Needles.leaves=='Needle-leaved'&obs4mod$Evergreen.Deciduous=='Evergreen']<-pchidx[4]
pch_site[obs4mod$Needles.leaves=='Broadleaved'&obs4mod$Evergreen.Deciduous=='Mixed']<-pchidx[5]
pch_site[obs4mod$Needles.leaves=='Mixed'&obs4mod$Evergreen.Deciduous=='Mixed']<-pchidx[6]

colidx<-c('blue','black','orange','red')
colleg<-c('boreal','temperate','mediterr.','tropical')
col_site<-numeric(nsite)
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


pdf('compare_sim_cruncep_vs_gfdb-plots.pdf')

plot_sitelocs(obs4mod$Longitude,obs4mod$Latitude,pch_site,col_site)

xvar<-obs4mod$BPo;yvar=model_mean$NPP*12.001;label='NPP (gC/m2/yr)'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$GPP;yvar=model_mean$GPP*12.001;label='GPP (gC/m2/yr)'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$BPEo;yvar=model_mean$NPP/model_mean$GPP;label='NPP/GPP' 
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$LAI;yvar=model_mean$LAI;label='LAI'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$Foliar.BM*0.48;yvar=model_mean$Leaf_C*12.001;label='Leaf C (gC/m2)'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$Fine.Root.BM*0.48;yvar=model_mean$FineRoot_C*12.001;label='Fine Root C (gC/m2)'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$Coarse.root.BM*0.48;yvar=model_mean$CoarseRoot_C*12.001;label='Coarse Root C (gC/m2)'
plot_mod2obs(xvar,yvar,label)

xvar<-obs4mod$Stem.BM
xvar[!is.na(obs4mod$Branch.BM)]<-xvar[!is.na(obs4mod$Branch.BM)]+obs4mod$Branch.BM[!is.na(obs4mod$Branch.BM)]
xvar[!is.na(obs4mod$Stump.BM)]<-xvar[!is.na(obs4mod$Stump.BM)]+obs4mod$Stump.BM[!is.na(obs4mod$Stump.BM)]
xvar<-xvar*0.48;yvar=(model_mean$SapWood_C+model_mean$HeartWood_C)*12.001;label='Woody C (gC/m2)'
plot_mod2obs(xvar,yvar,label)


dev.off()
