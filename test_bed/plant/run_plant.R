#!/usr/bin/env Rscript
# program to plot plant model output 
#  changed pattern to "*.txt" in line 8 [Jan Engel 20160315]
#  set pdf options [Jan Engel 20160316]

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

make=read.delim('../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=trim(make[trim(make[,1])=="COMP",2])
system(paste("../../",comp,"/bin/plant.x",sep=''))

#read timestep of output file
namelist=read.delim('namelist.plant',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="read_output_interval_pool",2])
flist_pool = Sys.glob(paste("pool*",tstep_pool,"*.txt",sep='')) 
flist_pool = c(flist_pool,Sys.glob(paste("veg_diagnostics*",tstep_pool,"*.txt",sep=''))) 
tstep_flux=trim(namelist[trim(namelist[,1])=="read_output_interval_flux",2])
flist_flux = Sys.glob(paste("flux*",tstep_flux,"*.txt",sep='')) 
flist=c(flist_pool,flist_flux)
nplot=length(flist)

pdf("plant_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
for (i in 1:nplot) {
  dat=read.table(flist[i],header=T)
  nvar=dim(dat)[2]
  mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.5,1,0,0))
  for (v in 1:nvar) plot(dat[,v],ylab=colnames(dat)[v],type='l')
}

dev.off()
