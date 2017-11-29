#!/usr/bin/env Rscript

# this script does: postprocessing of output data (e.g. plot plant model output)
# it is called either by contrib/site_level/run_slm.sh (HPC simulations) or run-slm_test-bed.sh

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])
soilLayer <- trim(namelist[trim(namelist[,1])=="nsoil_energy",2])
soilLayer <- as.numeric(soilLayer)

# message
#print(" postprocessing SLM plant output")

# get pool output files
flist_pool = Sys.glob(paste0("vegpool*","_",tstep_pool,"*.txt")) 
flist_pool = c(flist_pool,Sys.glob(paste0("soilpool*","_",tstep_pool,"*.txt"))) 
flist_pool = c(flist_pool,Sys.glob(paste0("veg_diagnostics*","_",tstep_pool,"*.txt")))
flist_pool = c(flist_pool,Sys.glob(paste0("soil_diagnostics*","_",tstep_pool,"*.txt")))
flist_pool = c(flist_pool,Sys.glob(paste0("ecopool*","_",tstep_pool,"*.txt")))  
# get flux output files  
flist_flux = Sys.glob(paste0("vegflux*","_",tstep_pool,"*.txt")) 
flist_flux = c(flist_flux,Sys.glob(paste0("soilflux*","_",tstep_pool,"*.txt"))) 

#create overall file list
flist=c(flist_pool,flist_flux)
nplot=length(flist)

# plot data
pdf("land_model_plots.pdf", paper="a4", height=10, width=8, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")

for (i in 1:nplot) {
  dat=read.table(flist[i],header=T,sep="")
  nvar=dim(dat)[2]
  if(!grepl("_lay",flist[i])) {
    print(paste('   reading',flist[i]))
    effColNum <- ncol(dat) -1 
    mat=matrix(seq(1,2*round(effColNum/2)+2,1),nrow=round(effColNum/2)+1,ncol=2,byrow=T) 
    layout(mat)
    par(mai=c(0.4,0.5,0.14,0),oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins
    xAxis <- tstep_pool
    if(grepl("flux",flist[i])) {
      xAxis <- tstep_flux
    }
    for (v in 2:nvar) {
      if(!all(is.na(dat[,v]))) plot(dat[,v],ylab=colnames(dat)[v],type='l',xlab=xAxis, tck=0.03, mgp=c(1.5, 0.2, 0))
    }
  }
}

garbage <- dev.off()


# plot data with information per soil layer
flistLayFiles <- grepl("_lay",flist[])

# if there is at least one TRUE in flistLayFiles
if (any(flistLayFiles)){

  pdf("land_model_plots_lay.pdf", paper="a4", height=10, width=8, pointsize=12, 
      bg="white", family="Helvetica", version="1.6")

  for (i in 1:nplot) {
    dat=read.table(flist[i],header=T,sep="")
    nvar=dim(dat)[2]
    if(grepl("_lay",flist[i])) {
      print(paste('   reading',flist[i]))
      rowNum <- nrow(dat)
      varCol <- ncol(dat) -1 
      varNum <- varCol/ soilLayer
      mat=matrix(seq(1,2*round(soilLayer/2)+2,1),nrow=round(soilLayer/2)+1,ncol=2,byrow=T) 
      bottomMargin <- 18
      if (soilLayer > 5) bottomMargin <- 6  # reduce bottom margin to leave enough space for plots if soil layer numers is large
      
      xAxis <- tstep_pool
      if(grepl("flux",flist[i])) {
        xAxis <- tstep_flux
      }
    
      countVar <- 0
      for (v in 1:varNum) {
        layout(mat)
        par(mai=c(0.6,0.8,0,0),oma=c(bottomMargin,1,1,1)) # oma = outer margins
        a <- soilLayer *countVar +2
        b <- a + soilLayer -1
        dat2 <- dat[,a:b]
        for (j in 1:soilLayer) {
          plot(1:rowNum, dat2[,j],ylab=colnames(dat2)[j],type='l',xlab=xAxis, mgp=c(1.84, 0.64, 0))
        }
        countVar <- countVar +1
      }
    }
  }
  garbage <- dev.off()
} # end: if (any(flistLayFiles))
