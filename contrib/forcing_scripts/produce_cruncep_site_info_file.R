#script to generate surface run files for SLM

#available met forcing
metsites<-read.table('siteinfo_noheader.dat',stringsAsFactors=F)
metsite_list<-metsites[,1]

# available site characteristics
dat<-read.csv('info_forest.csv',sep=';',stringsAsFactors=F)
pft1<-numeric(dim(dat)[1])
pft1[dat[,3]=='Boreal Humid']='B'
pft1[dat[,3]=='Boreal Semi-arid']='B'
pft1[dat[,3]=='Mediterranean Warm']='Te'
pft1[dat[,3]=='Temperate Humid']='Te'
pft1[dat[,3]=='Temperate Semi-arid']='Te'
pft1[dat[,3]=='Tropical Humid']='Tr'
pft1[dat[,3]=='Tropical Semi-arid']='Tr'

pft2<-numeric(dim(dat)[1])
pft2[dat[,4]=='Broadleaved']='B'
pft2[dat[,4]=='Mixed']='M'
pft2[dat[,4]=='Needle-leaved']='N'

pft3<-numeric(dim(dat)[1])
pft3[dat[,5]=='Deciduous']='S'
pft3[dat[,5]=='Mixed']='M'
pft3[dat[,5]=='Evergreen']='E'
pft3[pft1=='Tr'&pft3=='S']<-'R'

pft1[pft1=='Te'&pft2=='N'&pft3=='S']<-'B'
pft1[pft1=='Tr'&pft2=='N'&pft3=='E']<-'Te'


site<-dat[,1]
pft=paste(pft1,pft2,pft3,sep='')
lon_orig<-dat[,8]
lat_orig<-dat[,6]
sand<-dat[,14]/100
silt<-dat[,15]/100
awc<-dat[,13]
age<-dat[,10]
nstat<-dat[,11];nstat[is.na(nstat)]<-0

#get supplementary data from Soil Grids Database
library(ncdf)
ipath<-'/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_static/soilgrids/v0_5_1/Data/'
ncid<-open.ncdf(paste(ipath,'SNDPPT.soilgrid.3600.1800.7.nc',sep=''))
lon_sg<-get.var.ncdf(ncid,'longitude')
lat_sg<-get.var.ncdf(ncid,'latitude')
sand_sg<-get.var.ncdf(ncid,'SNDPPT')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'SLTPPT.soilgrid.3600.1800.7.nc',sep=''))
silt_sg<-get.var.ncdf(ncid,'SLTPPT')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'BLDFIE.soilgrid.3600.1800.7.nc',sep=''))
bulk_density_sg<-get.var.ncdf(ncid,'BLDFIE')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'PHIHOX.soilgrid.3600.1800.7.nc',sep=''))
ph_sg<-get.var.ncdf(ncid,'PHIHOX')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'AWCtS.soilgrid.3600.1800.7.nc',sep=''))
fc_vol_sg<-get.var.ncdf(ncid,'AWCtS')
close.ncdf(ncid)
ncid<-open.ncdf(paste(ipath,'WWP.soilgrid.3600.1800.7.nc',sep=''))
pwp_vol_sg<-get.var.ncdf(ncid,'WWP')
close.ncdf(ncid)
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
taxnwrb<-taxousda<-depth<-ph<-bulk_density<-numeric(58)
awc_sl<-array(NA,c(580,7))
for (s in 1:580) {
  lonidx<-which.min(abs(lon_sg-lon_orig[s]))
  latidx<-which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s]<-mean(bulk_density_sg[lonidx,latidx,])
  ph[s]<-mean(ph_sg[lonidx,latidx,])/10
  depth[s]<-depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]<-taxousda_sg[lonidx,latidx]
  taxnwrb[s]<-taxnwrb_sg[lonidx,latidx]

  if(is.na(sand[s])) sand[s]<-weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  if(is.na(silt[s])) silt[s]<-weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  awc_sl[s,]<-(fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  if(is.na(awc[s]))  awc[s]<-weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}

#read lithology map
ncid<-open.ncdf('/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc')
lon_glim<-get.var.ncdf(ncid,'lon')
lat_glim<-get.var.ncdf(ncid,'lat')
glim_class<-get.var.ncdf(ncid,'GLim')
glim_class_site<-numeric(580)
for (s in 1:580) {
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
age[is.na(age)]<-500
taxousda[is.na(taxousda)]<-30
taxnwrb[is.na(taxnwrb)]<-27
glim_class_site[is.na(glim_class_site)]<-2
plantyear<-2000-age

# available site information (GPP, NPP and sections of it, even if not matching)
gfdb_all<-read.csv('info_forest_DATA.csv',sep=';',stringsAsFactors=F)
all_idx<-!(is.na(gfdb_all$TNPP)&is.na(gfdb_all$GPP)&is.na(gfdb_all$ANPP))&gfdb_all[,1]!='MAN'&gfdb_all[,1]!='PUU'&gfdb_all[,1]!='SKO'&gfdb_all[,1]!='SUE'&gfdb_all[,1]!='YUS_01'&gfdb_all[,1]!='YUS_02'&gfdb_all[,1]!='CHAM_01'&gfdb_all[,1]!='CHAM_02'&gfdb_all[,1]!='CHAM_03'

dummy_selected<-gfdb_all[all_idx,]
dummy_selected_sitelist<-dummy_selected[,1]
not_mixed<-logical(length(dummy_selected_sitelist))
i=1
for (gfdb_site in dummy_selected_sitelist) {
  idx<-gfdb_site==site
  if(pft2[idx]!='M'&pft3[idx]!='M') {not_mixed[i]<-TRUE}
  else {not_mixed[i]<-FALSE}
  i<-i+1
}
gfdb_all_selected<-dummy_selected[not_mixed,]
gfdb_all_selected_sitelist<-gfdb_all_selected[,1]

#match GFDB data with metsite database and correct name convention
oldsite=""
count=1
metsite_list_found=c()
gfdb_sitename_new=numeric(length(gfdb_all_selected_sitelist))
i=1
for (gfdb_site in gfdb_all_selected_sitelist){
 dummy<-gfdb_site
 # renamed site list
 if(dummy=='BILY') dummy<-'BLY'
 if(dummy=='CASTEL') dummy<-'CAP'
 if(dummy=='CHAK') dummy<-'CHX'
 if(strtrim(dummy,4)=='CHAM') dummy<-'CHM'
 if(dummy=='FRAN') dummy<-'FRN'
 if(strtrim(dummy,4)=='HARD') dummy<-'HRD'
 if(dummy=='HEST') dummy<-'HST'
 if(dummy=='KOIN') dummy<-'KIN'
 if(dummy=='KOIW') dummy<-'KIW'
 if(dummy=='LAUPA') dummy<-'LPA'
 if(dummy=='LOOB') dummy<-'LOB'
 if(dummy=='PRI_SOAS') dummy<-'PRS'
 if(dummy=='PRI_SOBS') dummy<-'PRB'
 if(dummy=='PRI_SOJP') dummy<-'PRP'
 if(dummy=='PRI_SYJP') dummy<-'PRP'
 if(dummy=='THO_S') dummy<-'THS'
 if(dummy=='TOMA') dummy<-'TMA'
 if(dummy=='TURA') dummy<-'TRA'
 if(dummy=="THO_W13") {oldsite="THO";count=6} 
 if(dummy=='USA_01') dummy<-'US1'
 if(dummy=='USA_02') dummy<-'US2'
 if(dummy=='USA_03') dummy<-'US2'
 if(dummy=='USA_04') dummy<-'US2'
 if(dummy=='USA_05') dummy<-'US2'
 if(dummy=='WAL_01') dummy<-'WL0'
 if(dummy=='WAL_02') dummy<-'WL0'
 if(dummy=='WAL_03') dummy<-'WL0'
 if(dummy=='WAL_04') dummy<-'WL0'
 if(dummy=='WC_00') dummy<-'WCA'

 idx<-metsite_list==strtrim(dummy,3)
 if(any(idx)) {
    newsite<-metsite_list[idx]
print(c(dummy,newsite))
    if(newsite==oldsite) {count=count+1}
    else {count=1}
    if(count<10) {gfdb_sitename_new[i]<-paste(metsite_list[idx],'_0',count,sep='')}
    if(count>=10) {gfdb_sitename_new[i]<-paste(metsite_list[idx],'_',count,sep='')}
    oldsite=newsite
    i=i+1
 }
 else {print(dummy)}
}

i=1
siteinfo<-array(NA,c(length(gfdb_all_selected_sitelist),21))
sitedata<-array(NA,c(length(gfdb_all_selected_sitelist),length(gfdb_all_selected[1,])+1))
for (gfdb_site in gfdb_all_selected_sitelist) {
  idx<-gfdb_site==site
  idx2<-gfdb_site==gfdb_all_selected_sitelist
    siteinfo[i,1]<-gfdb_sitename_new[idx2]
    siteinfo[i,2]<-lon_orig[idx]
    siteinfo[i,3]<-lat_orig[idx]
    siteinfo[i,4]<-pft[idx]
    siteinfo[i,5]<-1901
    siteinfo[i,6]<-2015
    siteinfo[i,7]<-lon_orig[idx]
    siteinfo[i,8]<-clay[idx]
    siteinfo[i,9]<-silt[idx]
    siteinfo[i,10]<-sand[idx]
    siteinfo[i,11]<-awc[idx]
    siteinfo[i,12]<-bulk_density[idx]
    siteinfo[i,13]<-ph[idx]
    siteinfo[i,14]<-taxousda[idx]
    siteinfo[i,15]<-taxnwrb[idx]
    siteinfo[i,16]<-glim_class_site[idx]
    siteinfo[i,17]<-(-9999) #LAI
    siteinfo[i,18]<-(-9999) #Nleaf
    siteinfo[i,19]<-(-9999) #SLA
    siteinfo[i,20]<-(-9999) #Height
    siteinfo[i,21]<-plantyear[idx]

    sitedata[i,1]<-gfdb_sitename_new[idx2]
    for (j in 2:11) sitedata[i,j]<-gfdb_all_selected[idx2,j]
    sitedata[i,12]<-gfdb_all_selected[idx2,1]
    i<-i+1
}

colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end','lon_gmt','clay','silt','sand','awc','bd','ph','taxusda','taxnwrb','lith_glim','LAI','Nleaf','SLA','Height','PlantYear')
write.table(siteinfo,file='all_gfdb_sites_1901-2015_list.dat',quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
write.table(siteinfo,file='all_gfdb_sites_1984-2013_list.dat',quote=F,row.names=F,col.names=T)
colnames(sitedata)<-c(colnames(gfdb_all_selected),'Original.Site.Code')
write.csv(sitedata,file='global_forest_database_all_sites_01082017.csv',quote=F,row.names=F,col.names=T)

# available site information (CUE and compontents
gfdb_cue<-read.csv('global_forest_database_21092016.csv',sep=',',stringsAsFactors=F)
all_idx<-gfdb_cue[,1]!='XIS_01'&gfdb_cue[,1]!='XIS_02'
dummy_selected<-gfdb_cue[all_idx,]
dummy_selected_sitelist<-dummy_selected[,1]
not_mixed<-logical(length(dummy_selected_sitelist))
i=1
for (gfdb_site in dummy_selected_sitelist) {
  idx<-gfdb_site==site
  if(gfdb_site=='CAL') idx<-'REN'==site
  if(gfdb_site=='NAU_M') idx<-'NAU'==site
  if(gfdb_site=='NAU_U') idx<-'NAU'==site
  if(gfdb_site=='POPF') idx<-'POP_01'==site
  if(gfdb_site=='THO_D15') idx<-'THO_D12'==site
  if(gfdb_site=='THO_D71') idx<-'THO_71'==site
  if(pft2[idx]!='M'&pft3[idx]!='M') {not_mixed[i]<-TRUE}
  else {not_mixed[i]<-FALSE}
  i<-i+1
}
gfdb_cue<-dummy_selected[not_mixed,]
#match GFDB data with metsite database and correct name convention
oldsite=""
count=1
gfdb_cue_sitelist<-gfdb_cue[,1]
gfdb_cue_sitename_new=numeric(length(gfdb_cue_sitelist))
i=1
for (gfdb_site in gfdb_cue_sitelist){
 dummy<-gfdb_site
 # renamed site list
 if(dummy=='BILY') dummy<-'BLY'
 if(dummy=='CAL') dummy<-'REN' #only for met purposes!!
 if(dummy=='CASTEL') dummy<-'CAP'
 if(dummy=='CHAK') dummy<-'CHX'
 if(strtrim(dummy,4)=='CHAM') dummy<-'CHM'
 if(dummy=='FRAN') dummy<-'FRN'
 if(strtrim(dummy,4)=='HARD') dummy<-'HRD'
 if(dummy=='HEST') dummy<-'HST'
 if(dummy=='KOIN') dummy<-'KIN'
 if(dummy=='KOIW') dummy<-'KIW'
 if(dummy=='LAUPA') dummy<-'LPA'
 if(dummy=='LOOB') dummy<-'LOB'
 if(dummy=='PRI_SOAS') dummy<-'PRS'
 if(dummy=='PRI_SOBS') dummy<-'PRB'
 if(dummy=='PRI_SOJP') dummy<-'PRP'
 if(dummy=='PRI_SYJP') dummy<-'PRP'
 if(dummy=='THO_S') dummy<-'THS'
 if(dummy=='TOMA') dummy<-'TMA'
 if(dummy=='TURA') dummy<-'TRA'
 if(dummy=="THO_W13") {oldsite="THO";count=7} 
 if(dummy=='USA_01') dummy<-'US1'
 if(dummy=='USA_02') dummy<-'US2'
 if(dummy=='USA_03') dummy<-'US2'
 if(dummy=='USA_04') dummy<-'US2'
 if(dummy=='USA_05') dummy<-'US2'
 if(dummy=='WAL_01') dummy<-'WL0'
 if(dummy=='WAL_02') dummy<-'WL0'
 if(dummy=='WAL_03') dummy<-'WL0'
 if(dummy=='WAL_04') dummy<-'WL0'
 if(dummy=='WC_00') dummy<-'WCA'

 idx<-metsite_list==strtrim(dummy,3)
 if(any(idx)) {
    newsite<-metsite_list[idx]
    if(newsite==oldsite) {count=count+1}
    else {count=1}
    if(count<10) {gfdb_cue_sitename_new[i]<-paste(metsite_list[idx],'_0',count,sep='')}
    if(count>=10) {gfdb_cue_sitename_new[i]<-paste(metsite_list[idx],'_',count,sep='')}
    oldsite=newsite
    i=i+1
  } 
 else {print(dummy)}
}

i=1
siteinfo<-array(NA,c(length(gfdb_cue_sitelist),21))
sitedata<-array(NA,c(length(gfdb_cue_sitelist),length(gfdb_cue[1,])+1))
for (gfdb_site in gfdb_cue_sitelist) {
  idx<-gfdb_site==site
  if(gfdb_site=='CAL') idx<-'REN'==site
  if(gfdb_site=='NAU_M') idx<-'NAU'==site
  if(gfdb_site=='NAU_U') idx<-'NAU'==site
  if(gfdb_site=='POPF') idx<-'POP_01'==site
  if(gfdb_site=='THO_D15') idx<-'THO_D12'==site
  if(gfdb_site=='THO_D71') idx<-'THO_71'==site
  idx2<-gfdb_site==gfdb_cue_sitelist
  siteinfo[i,1]<-gfdb_cue_sitename_new[idx2]
    siteinfo[i,2]<-lon_orig[idx]
    siteinfo[i,3]<-lat_orig[idx]
    siteinfo[i,4]<-pft[idx]
    siteinfo[i,5]<-1901
    siteinfo[i,6]<-2015
    siteinfo[i,7]<-lon_orig[idx]
    siteinfo[i,8]<-clay[idx]
    siteinfo[i,9]<-silt[idx]
    siteinfo[i,10]<-sand[idx]
    siteinfo[i,11]<-awc[idx]
    siteinfo[i,12]<-bulk_density[idx]
    siteinfo[i,13]<-ph[idx]
    siteinfo[i,14]<-taxousda[idx]
    siteinfo[i,15]<-taxnwrb[idx]
    siteinfo[i,16]<-glim_class_site[idx]
    siteinfo[i,17]<-(-9999) #LAI
    siteinfo[i,18]<-(-9999) #Nleaf
    siteinfo[i,19]<-(-9999) #SLA
    siteinfo[i,20]<-(-9999) #Height
    siteinfo[i,21]<-plantyear[idx]

  sitedata[i,1]<-gfdb_cue_sitename_new[idx2]
  for (j in 2:length(gfdb_cue[1,])) sitedata[i,j]<-gfdb_cue[idx2,j]
  sitedata[i,length(gfdb_cue[1,])+1]<-gfdb_cue[idx2,1]
  i<-i+1
}
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end','lon_gmt','clay','silt','sand','awc','bd','ph','taxusda','taxnwrb','lith_glim','LAI','Nleaf','SLA','Height','PlantYear')
write.table(siteinfo,file='cue_gfdb_sites_1901-2015_list.dat',quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
write.table(siteinfo,file='cue_gfdb_sites_1984-2013_list.dat',quote=F,row.names=F,col.names=T)
colnames(sitedata)<-c(colnames(gfdb_cue),'Original.Site.Code')
write.csv(sitedata,file='global_forest_database_cue_sites_01082017.csv',quote=F,row.names=F,col.names=T)

#SPP sites
spp_site_list<-c('BBR','LUE')
lon_orig<-lat_orig<-numeric(2)
i=1
for (spp_site in spp_site_list) {
  idx<-spp_site==metsite_list
  lon_orig[i]<-metsites[idx,2]
  lat_orig[i]<-metsites[idx,3]
  i<-i+1
}

awc<-clay<-silt<-sand<-taxnwrb<-taxousda<-depth<-ph<-bulk_density<-numeric(2)
awc_sl<-array(NA,c(2,7))
for (s in 1:2) {
  lonidx<-which.min(abs(lon_sg-lon_orig[s]))
  latidx<-which.min(abs(lat_sg-lat_orig[s]))
  bulk_density[s]<-mean(bulk_density_sg[lonidx,latidx,])
  #ph[s]<-mean(ph_sg[lonidx,latidx,])/10
  depth[s]<-depth_to_bedrock_sg[lonidx,latidx]/100
  taxousda[s]<-taxousda_sg[lonidx,latidx]
  taxnwrb[s]<-taxnwrb_sg[lonidx,latidx]

  #sand[s]<-weighted.mean(sand_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  #silt[s]<-weighted.mean(silt_sg[lonidx,latidx,],depth_sl,na.rm=TRUE)/100
  awc_sl[s,]<-(fc_vol_sg[lonidx,latidx,]-pwp_vol_sg[lonidx,latidx,])/100
  awc[s]<-weighted.mean(awc_sl[s,],depth_sl,na.rm=TRUE)*depth[s]*1000
}
sand<-c(22.59,76.6)/100
silt<-c(51.31,17.13)/100
ph<-c(5.2,4.3)
clay<-1-sand-silt

#read lithology map
ncid<-open.ncdf('/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_static/GLiM/v1_0/Data/GLim.720.360.nc')
lon_glim<-get.var.ncdf(ncid,'lon')
lat_glim<-get.var.ncdf(ncid,'lat')
glim_class<-get.var.ncdf(ncid,'GLim')
glim_class_site<-numeric(2)
for (s in 1:2) {
  lonidx<-which.min(abs(lon_glim-lon_orig[s]))
  latidx<-which.min(abs(lat_glim-lat_orig[s]))
  glim_class_site[s]<-glim_class[lonidx,latidx]
}

#replace NAN iwth default values
siteinfo<-array(NA,c(2,21))
i=1
for (spp_site in spp_site_list) {
  idx<-spp_site==metsite_list
  siteinfo[i,1]<-paste(spp_site,'_01',sep='')
    siteinfo[i,2]<-metsites[idx,2]
    siteinfo[i,3]<-metsites[idx,3]
    siteinfo[i,4]<-'TeBS'
    siteinfo[i,5]<-1901
    siteinfo[i,6]<-2015
    siteinfo[i,7]<-metsites[idx,2]
    siteinfo[i,8]<-clay[i]
    siteinfo[i,9]<-silt[i]
    siteinfo[i,10]<-sand[i]
    siteinfo[i,11]<-awc[i]
    siteinfo[i,12]<-bulk_density[i]
    siteinfo[i,13]<-ph[i]
    siteinfo[i,14]<-taxousda[i]
    siteinfo[i,15]<-taxnwrb[i]
    siteinfo[i,16]<-glim_class_site[i]
    siteinfo[i,17]<-(-9999) #LAI
    siteinfo[i,18]<-(-9999) #Nleaf
    siteinfo[i,19]<-(-9999) #SLA
    siteinfo[i,20]<-(-9999) #Height
    siteinfo[i,21]<-1890
  i=i+1
}
colnames(siteinfo)<-c('Site-ID','lon','lat','PFT','start','end','lon_gmt','clay','silt','sand','awc','bd','ph','taxusda','taxnwrb','lith_glim','LAI','Nleaf','SLA','Height','PlantYear')
write.table(siteinfo,file='spp_sites_1901-2015_list.dat',quote=F,row.names=F,col.names=T)
siteinfo[,5]<-1984
siteinfo[,6]<-2013
write.table(siteinfo,file='spp_sites_1984-2013_list.dat',quote=F,row.names=F,col.names=T)



