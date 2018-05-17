#script to generate surface run files for SLM
library(ncdf)

ipath1='/Net/Groups/BGI/model/OCN/input/point/FACE/'
opath="/Net/Groups/BGI/data/Model/SLM1/input/point/FACE/"

allsiteinfo<-read.table(paste('site_list.dat',sep=''),stringsAsFactor=F,header=T)
allsitestartyear<-allsiteinfo[,2]
allsiteendyear<-allsiteinfo[,3]
allsitename<-allsiteinfo[,1]
allsite_newname<-allsiteinfo[,6]
site_list<-allsitename
pft<-allsiteinfo[,4]

lon_orig<-lat_orig<-numeric(length(site_list))
for (s in 1:length(site_list)) {
  fid<-open.ncdf(paste(ipath1,site_list[s],'_FACE.nc',sep=''))
    lon_orig[s]<-get.var.ncdf(fid,"x")
    lat_orig[s]<-get.var.ncdf(fid,"y")
}

gmt_ref<-allsiteinfo[,15]
sand<-allsiteinfo[,12]
silt<-allsiteinfo[,11]
clay<-1-sand-silt
awc<-allsiteinfo[,13]
ph<-allsiteinfo[,14]
lai<-allsiteinfo[,10]*0+-9999
nleaf<-allsiteinfo[,11]*0+-9999
sla<-allsiteinfo[,12]*0+-9999
height<-allsiteinfo[,13]*0+-9999
plantyear<-allsiteinfo[,16]

#get supplementary data from Soil Grids Database
ipath<-'/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/'
ncid<-open.ncdf(paste(ipath,'SNDPPT.soilgrid.3600.1800.7.nc',sep=''))
lon_sg<-get.var.ncdf(ncid,'longitude')
lat_sg<-get.var.ncdf(ncid,'latitude')
#sand_sg<-get.var.ncdf(ncid,'SNDPPT')
#close.ncdf(ncid)
#ncid<-open.ncdf(paste(ipath,'SLTPPT.soilgrid.3600.1800.7.nc',sep=''))
#silt_sg<-get.var.ncdf(ncid,'SLTPPT')
#close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'BLDFIE.soilgrid.3600.1800.7.nc',sep=''))
bulk_density_sg<-get.var.ncdf(ncid,'BLDFIE')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'PHIHOX.soilgrid.3600.1800.7.nc',sep=''))
ph_sg<-get.var.ncdf(ncid,'PHIHOX')
close.ncdf(ncid)
#ncid<-open.ncdf(paste(ipath,'AWCtS.soilgrid.3600.1800.7.nc',sep=''))
#fc_vol_sg<-get.var.ncdf(ncid,'AWCtS')
#close.ncdf(ncid)
#ncid<-open.ncdf(paste(ipath,'WWP.soilgrid.3600.1800.7.nc',sep=''))
#pwp_vol_sg<-get.var.ncdf(ncid,'WWP')
#close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'BDRICM.soilgrid.3600.1800.nc',sep=''))
depth_to_bedrock_sg<-get.var.ncdf(ncid,'BDRICM')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'TAXOUSDA_10km.soilgrid.3600.1800.nc',sep=''))
taxousda_sg<-get.var.ncdf(ncid,'TAXOUSDA_10km')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'TAXNWRB.soilgrid.3600.1800.nc',sep=''))
taxnwrb_sg<-get.var.ncdf(ncid,'TAXNWRB')
close.ncdf(ncid)

depth_begin<-c(0.0,0.05,0.15,0.3,0.6,1.0,2.0)
depth_end<-c(0.05,0.15,0.3,0.6,1.0,2.0,100)
depth_sl<-depth_end-depth_begin
taxnwrb<-taxousda<-depth<-bulk_density<-numeric(length(site_list))
awc_sl<-array(NA,c(length(site_list),7))
for (s in 1:length(site_list)) {
  lonidx<-which.min(abs(lon_sg-lon_orig[s]))
  latidx<-which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s]<-mean(bulk_density_sg[lonidx,latidx,])
  if(ph[s]==6.5) ph[s]<-mean(ph_sg[lonidx,latidx,])/10
#  depth[s]<-depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]<-taxousda_sg[lonidx,latidx]
  taxnwrb[s]<-taxnwrb_sg[lonidx,latidx]

#  if(is.na(sand[s])) sand[s]<-weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
#  if(is.na(silt[s])) silt[s]<-weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
#  awc_sl[s,]<-(fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
#  if(is.na(awc[s]))  awc[s]<-weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}

#read lithology map
ncid<-open.ncdf('/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc')
lon_glim<-get.var.ncdf(ncid,'lon')
lat_glim<-get.var.ncdf(ncid,'lat')
glim_class<-get.var.ncdf(ncid,'GLim')
glim_class_site<-numeric(length(site_list))
for (s in 1:length(site_list)) {
  lonidx<-which.min(abs(lon_glim-lon_orig[s]))
  latidx<-which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s]<-glim_class[lonidx,latidx]
}

#replace NAN iwth default values
sand[is.na(sand)]<-0.4
silt[is.na(silt)]<-0.4
clay<-1-sand-silt
bulk_density[is.na(bulk_density)]<-1500
ph[is.na(ph)]<-6
awc[is.na(awc)]<-200
taxousda[is.na(taxousda)]<-30
taxnwrb[is.na(taxnwrb)]<-27
glim_class_site[is.na(glim_class_site)]<-2

siteinfo<-array(NA,c(length(site_list),21))
siteinfo[,1]<-allsite_newname
siteinfo[,2]<-lon_orig
siteinfo[,3]<-lat_orig
siteinfo[,4]<-pft
siteinfo[,5]<-allsitestartyear
siteinfo[,6]<-allsiteendyear
siteinfo[,7]<-gmt_ref
siteinfo[,8]<-clay
siteinfo[,9]<-silt
siteinfo[,10]<-sand
siteinfo[,11]<-awc
siteinfo[,12]<-bulk_density
siteinfo[,13]<-ph
siteinfo[,14]<-taxousda
siteinfo[,15]<-taxnwrb
siteinfo[,16]<-glim_class_site
siteinfo[,17]<-lai
siteinfo[,18]<-nleaf
siteinfo[,19]<-sla
siteinfo[,20]<-height
siteinfo[,21]<-plantyear

colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end','lon_gmt','clay','silt','sand','awc','bd','ph','taxusda','taxnwrb','lith_glim','LAI','Nleaf','SLA','Height','PlantYear')
write.table(siteinfo,file='all_sites_1996-2008_list.dat',quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1901
siteinfo[,6]<-2015
write.table(siteinfo,file='all_sites_1901-2015_list.dat',quote=F,row.names=F,col.names=T)

