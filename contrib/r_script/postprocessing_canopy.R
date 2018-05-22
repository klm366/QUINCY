#!/usr/bin/env Rscript

# postprocessing of output data (e.g. plot canopy model output)
# called from contrib/site_level/run_slm.sh (HPC simulations) or run_canopy.R (local SLM runs)

###
### copy pasted from the postprocessing_plant.R to cover the flux output files
###
print("postprocessing SLM output")

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# message
#print(" postprocessing SLM plant output")

# get flux output files  
flist_flux = Sys.glob(paste0("vegflux*","_",tstep_flux,"*.txt")) 
# get other output files
flist_other = Sys.glob(paste0("atm_forcing","_",tstep_flux,"*.txt"))

#create overall file list
flist=c(flist_flux,flist_other)
nplot=length(flist)

# plot data
pdf("canopy_model_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")

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

garbage <- dev.off()



# message
#print(" finished postprocessing SLM output !")
