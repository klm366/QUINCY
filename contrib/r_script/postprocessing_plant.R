#!/usr/bin/env Rscript

# postprocessing of output data from QUINCY plant model
# it is called either by contrib/site_level/run_slm.sh (HPC simulations) or contrib/test_bed/run-slm_test-bed.sh

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# get pool output files
flist_pool = Sys.glob(paste0("vegpool*","_",tstep_pool,"*.txt")) 
flist_pool = c(flist_pool,Sys.glob(paste0("soilpool*","_",tstep_pool,"*.txt"))) 
flist_pool = c(flist_pool,Sys.glob(paste0("veg_diagnostics*","_",tstep_pool,"*.txt"))) 
flist_pool = c(flist_pool,Sys.glob(paste0("ecopool*","_",tstep_pool,"*.txt")))  
# get flux output files  
flist_flux = Sys.glob(paste0("vegflux*","_",tstep_pool,"*.txt"))
# get other output files
flist_other = Sys.glob(paste0("atm_forcing","_",tstep_flux,"*.txt"))

#create overall file list
flist=c(flist_pool,flist_flux,flist_other)
nplot=length(flist)

# plot data
pdf("plant_model_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")

for (i in 1:nplot) {
  print(paste('   reading',flist[i]))
  # only output files that do not contain layer information
  if(!grepl("_lay",flist[i])) {
    dat=read.table(flist[i],header=T)
    nvar=dim(dat)[2]
    mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
    layout(mat)
    par(mai=c(0.4,0.5,0.14,0),oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins
    xAxis <- tstep_pool
    if(grepl("flux",flist[i])) {
      xAxis <- tstep_flux
    }
    for (v in 2:nvar) {
      if(!all(is.na(dat[,v]))) plot(dat[,v], ylab=colnames(dat)[v], type='l', xlab=xAxis, tck=0.02, mgp=c(1.5, 0.2, 0))
    }
  }
}

garbage <- dev.off()
