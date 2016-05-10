#!/usr/bin/env Rscript
# program to plot plant model output, for multiple tiles separately (if applicable)

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# message
print(" running model simulations !")

# run plant.x
make=read.delim('../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=trim(make[trim(make[,1])=="COMP",2])
system(paste("../../",comp,"/bin/plant.x",sep=''))

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.plant',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# get number of tiles calculated -> needs "poolC_tilex_weekly.txt" as output files
# list all files in folder, grep for "poolC_tile*", and check vector length
tileNum<-length(grep("poolC_tile*",list.files()))

# message
print(" processing output data")

# loop over each tiles:
for(ti in 1:tileNum)
{
  print(paste0(" tile ",ti,"/",tileNum)) # message
  flist_pool = Sys.glob(paste0("pool*_tile",ti,"_",tstep_pool,".txt")) 
  flist_pool = c(flist_pool,Sys.glob(paste0("veg_diagnostics*_tile",ti,"_",tstep_pool,".txt"))) 

  flist_flux = Sys.glob(paste0("flux*_tile",ti,"_",tstep_pool,".txt")) 
  #create overall file list
  flist=c(flist_pool,flist_flux)
  nplot=length(flist)

  # plot data
  pdf(paste0("plant_plots_tile",ti,".pdf"), paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
  for (i in 1:nplot) {
    dat=read.table(flist[i],header=T)
    nvar=dim(dat)[2]
    mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
    layout(mat)
    par(mai=c(0.5,1,0,0))
    for (v in 1:nvar) plot(dat[,v],ylab=colnames(dat)[v],type='l')
  }

dev.off()
} # end loop over tiles

# message
print(" finished !")
