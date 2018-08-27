#!/usr/bin/env Rscript

## 
## create one .pdf file with figures from one specific simulation-output file for all sites
## 
## simply run this script within your simulation experiment folder "exp_*" to get another overview file
##   Rscript --vanilla ./postprocessing_site_overview.R
## 

# define the output file name without "_[timestep].txt"  (currently not working for _lay files with data per soil layer)
# e.g. vegpoolC  veg_diagnostics  soilfluxdN  ...
simOutFile <- "vegpoolC"


###################################################################################################

print(" Create one .pdf file with figs from one output file for all sites; default is vegpoolC .")

## detect sites within the folder this script is located
# get all directories names
dirnames  <- system('ls -d */',intern=T)
# remove directories from list that contain other stuff than the simulation results
dirnames  <- dirnames[!grepl("x_binary",dirnames)]
dirnames  <- dirnames[!grepl("parameter_list",dirnames)]
dirnames  <- dirnames[!grepl("default_files",dirnames)]
dirnames  <- dirnames[!grepl("parameter_sensitivity_files",dirnames)]
dirnames  <- dirnames[!grepl("out",dirnames)]
dirnames  <- dirnames[!grepl("error",dirnames)]
# number of directories/sites
nsite     <- length(dirnames)


## get timestep string from namelist of 1st site in list
namelistPath <- paste0(dirnames[1], "input/namelist.slm")
# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# read timestep of output file from namelist for pool & flux
namelist   <- read.delim(namelistPath, header=F, sep="=", stringsAsFactors=F)
tstep_pool <- trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux <- trim(namelist[trim(namelist[,1])=="output_interval_flux",2])

# set pool/flux timestep and x-axis label
tstep      <- tstep_pool
xAxis      <- tstep_pool
if(grepl("flux", simOutFile)) {
  tstep <- tstep_flux
  xAxis <- tstep_flux
}

## open output file
pdf(paste0(simOutFile, "_site_overview.pdf"), paper="a4", height=10, width=8, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")

## loop over sites
for (s in 1:nsite) {
  print(paste(s,dirnames[[s]]))
  site_str <- tail(unlist(strsplit(dirnames[[s]],'/')),n=1)
  
  dat <- read.table(paste0(dirnames[[s]], "output/", simOutFile, "_", tstep, ".txt"), header=T)

  # number of variables
  nvar=dim(dat)[2]
  # number of columns with variables, excluding the 1st row "Time"
  effColNum <- ncol(dat) -1 
  mat=matrix(seq(1,2*round(effColNum/2)+2,1),nrow=round(effColNum/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.4,0.5,0.14,0),oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins

  for (v in 2:nvar) {
    if(v < 4) maintext <- site_str
    else maintext <- ""
    if(!all(is.na(dat[,v]))){
      plot(dat[,v],ylab=colnames(dat)[v],type='l',xlab=xAxis, tck=0.03, mgp=c(1.5, 0.2, 0), main=maintext)
      # plot a brown horizontal line at 0 for sb_massBalance only
      if(grepl("sb_massBalance", simOutFile)) abline(h=0,col="sienna",lty=1,lwd=0.4)
    }
  }
}


garbage <- dev.off()
