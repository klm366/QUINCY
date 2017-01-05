#!/usr/bin/env Rscript

# this script does: postprocessing of output data (e.g. plot plant model output)
# it is called either by contrib/site_level/run_slm.sh (HPC simulations) or run_plant.R (local SLM runs)

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# get number of tiles calculated -> needs "poolC_tilex_weekly.txt" as output files
# list all files in folder, grep for "poolC_tile*", and check vector length
tileNum<-length(grep("vegpoolC_tile*",list.files()))
# if output files with "timestep" intervall exist, assume these are for "special case: fluxnet output" and do not consider for figure creation
filenameWithAnyNumber<-c("vegpoolC_tile","[0-9]+","_timestep.txt")
if (length( grep(paste0(filenameWithAnyNumber,collapse = ""),list.files()) ) > 0)
{
  tileNum <- tileNum - (tileNum/2)
}

# message
#print(" postprocessing SLM output")

# loop over each tiles:
for(ti in 1:tileNum)
{
  # message
  print(paste0(" tile ",ti,"/",tileNum))
  flist_pool = Sys.glob(paste0("vegpool*_tile",ti,"_",tstep_pool,"*.txt")) 
  flist_pool = c(flist_pool,Sys.glob(paste0("soilpool*_tile",ti,"_",tstep_pool,"*.txt"))) 
  flist_pool = c(flist_pool,Sys.glob(paste0("veg_diagnostics*_tile",ti,"_",tstep_pool,"*.txt"))) 

  flist_flux = Sys.glob(paste0("vegflux*_tile",ti,"_",tstep_pool,"*.txt")) 
  #create overall file list
  flist=c(flist_pool,flist_flux)
  nplot=length(flist)

  # plot data
  pdf(paste0("plant_plots_tile",ti,".pdf"), paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
  for (i in 1:nplot) {
    print(paste('   reading',flist[i]))
    dat=read.table(flist[i],header=T)
    nvar=dim(dat)[2]
    mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
    layout(mat)
    par(mai=c(0.5,1,0,0))
    for (v in 2:nvar) {
      if(!all(is.na(dat[,v]))) plot(dat[,v],ylab=colnames(dat)[v],type='l')
   }
  }

dev.off()
} # end loop over tiles

# message
#print(" finished postprocessing SLM output !")
