#!/usr/bin/env Rscript

##    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
##    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
##    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
##
##    This file is part of the QUINCY model.
##    The QUINCY model is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License version 3,
##    as published by the Free Software Foundation.
##
##    The QUINCY model is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.

# postprocessing of output data from QUINCY canopy model
# it is called either by contrib/site_level/run_slm.sh (HPC simulations) or contrib/test_bed/run-slm_test-bed.sh

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F,skip=17)
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# get flux output files  
flist_flux = Sys.glob(paste0("vegflux*","_",tstep_flux,"*.txt")) 
# get other output files
flist_other = Sys.glob(paste0("atm_forcing","_",tstep_flux,"*.txt"))

#create overall file list
flist  <- c(flist_flux,flist_other)
# remove output files that contain data from spinup
spinup <- which(!grepl("spinup",flist))
flist  <- flist[spinup]
# number of output files
nplot=length(flist)

# plot data
pdf("canopy_model_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")

for (i in 1:nplot) {
  print(paste('   reading ',flist[i]))
  dat=read.table(flist[i],header=T)
  nvar=dim(dat)[2]
  mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.4,0.5,0.14,0), oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins
  xAxis <- tstep_flux
  for (v in 2:nvar) {
    if(!all(is.na(dat[,v]))) plot(dat[,v], ylab=colnames(dat)[v], type='l', xlab=xAxis, tck=0.02, mgp=c(1.5, 0.2, 0))
  }
}

garbage <- dev.off()
