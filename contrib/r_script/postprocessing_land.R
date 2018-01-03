#!/usr/bin/env Rscript

# this script does: postprocessing of output data (e.g. plot plant model output)
# it is called either by contrib/site_level/run_slm.sh (HPC simulations) or run-slm_test-bed.sh

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# read timestep of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])
soilLayer <- trim(namelist[trim(namelist[,1])=="nsoil_energy",2])
soilLayer <- as.numeric(soilLayer)
soilLSeq  <- seq(1,soilLayer,1)

# get pool output files
flist_pool = Sys.glob(paste0("vegpool*","_",tstep_pool,"*.txt")) 
flist_pool = c(flist_pool,Sys.glob(paste0("soilpoolPh*","_",tstep_pool,"*.txt"))) 
flist_pool = c(flist_pool,Sys.glob(paste0("soilpool*","_",tstep_pool,"*.txt"))) 
flist_pool = c(flist_pool,Sys.glob(paste0("veg_diagnostics*","_",tstep_pool,"*.txt")))
flist_pool = c(flist_pool,Sys.glob(paste0("soil_diagnostics*","_",tstep_pool,"*.txt")))
flist_pool = c(flist_pool,Sys.glob(paste0("ecopool*","_",tstep_pool,"*.txt")))  
flist_pool = unique(flist_pool)
# get flux output files  
flist_flux = Sys.glob(paste0("vegflux*","_",tstep_flux,"*.txt"))
flist_flux = c(flist_flux,Sys.glob(paste0("soilflux*","_",tstep_flux,"*.txt")))
flist_flux = c(flist_flux,Sys.glob(paste0("mycoflux*","_",tstep_flux,"*.txt")))

# create overall file list
flist=c(flist_pool,flist_flux)
# number of output files
nOutFiles=length(flist)


###############################################################################
# plot data that is not "per soil layer data"
###############################################################################
pdf("land_model_plots.pdf", paper="a4", height=10, width=8, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")
print("Creating: land_model_plots.pdf")
# loop over all output files
for (i in 1:nOutFiles) {
  # only output files that do not contain layer information
  if(!grepl("_lay",flist[i])) {
    print(paste('   reading',flist[i]))
    dat=read.table(flist[i],header=T,sep="")
    nvar=dim(dat)[2]
    # number of columns with variables, excluding the 1st row "Time"
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


###############################################################################
# plot data with information per soil layer
###############################################################################
flistLayFiles <- grepl("_lay",flist[])

# if there is at least one TRUE in flistLayFiles
if (any(flistLayFiles)){

  pdf("land_model_plots_lay.pdf", paper="a4", height=10, width=8, pointsize=12, 
      bg="white", family="Helvetica", version="1.6")
  print("Creating: land_model_plots_lay.pdf")
  # loop over all output files
  for (i in 1:nOutFiles) {
    # only output files that contain layer information
    if(grepl("_lay",flist[i])) {
      print(paste('   reading',flist[i]))
      dat=read.table(flist[i],header=T,sep="")
      nvar=dim(dat)[2]
      # number of rows with data
      rowNum <- nrow(dat)
      # number of columns with variables, excluding the 1st row "Time"
      varCol <- ncol(dat) -1 
      # number of variables
      varNum <- varCol/ soilLayer
      mat=matrix(seq(1,2*round(soilLayer/2)+2,1),nrow=round(soilLayer/2)+1,ncol=2,byrow=T) 
      bottomMargin <- 18
      if (soilLayer > 8) bottomMargin <- 8     # reduce bottom margin to leave enough space for plots if soil layer numers is large
      if (soilLayer > 14) bottomMargin <- 0.1  # reduce bottom margin to leave enough space for plots if soil layer numers is very large
      
      xAxis <- tstep_pool
      if(grepl("flux",flist[i])) {
        xAxis <- tstep_flux
      }
    
      countVar <- 0
      for (v in 1:varNum) {
        layout(mat)
        par(mai=c(0.3,0.4,0,0),oma=c(bottomMargin,0.2,1.2,0.2)) # oma = outer margins
        a <- soilLayer *countVar +2
        b <- a + soilLayer -1
        dat2 <- dat[,a:b]
        # variable name
        vName<-colnames(dat2[1])
        vName<-substr(vName, 1, nchar(vName)-3)
        for (j in soilLSeq) {
          if (j < 10){ yAxisLabel<-substr(colnames(dat2)[j], nchar(colnames(dat2)[j])-1, nchar(colnames(dat2)[j]))
          } else       yAxisLabel<-substr(colnames(dat2)[j], nchar(colnames(dat2)[j])-2, nchar(colnames(dat2)[j]))
          plot(1:rowNum, dat2[,j], 
            ylab=yAxisLabel, xlab=xAxis, 
            type='l',
            tck=0.02, mgp=c(1.2, 0.14, 0))
          if (j == 1) mtext(vName, cex=0.9, side=3, line=0.04, adj=0,col="red3")
          if (j == 2) mtext(flist[i], cex=0.8, side=3, line=0.04, adj=0)
        }
        countVar <- countVar +1
      }
    }
  }
  garbage <- dev.off()
} # end: if (any(flistLayFiles))


###############################################################################
# plot soil profiles (all var in *_lay* output files)
###############################################################################
# if there is at least one TRUE in flistLayFiles
if (any(flistLayFiles)){

  # row numbers that cover the last 10 years, or all if simulation time less than 10 years
  # possible output periods: timestep daily weekly monthly yearly
  # not differentiating pool vs flux tstep_* !
  timeFac <- 0    # scaling factor to get the correct number of rows
  spanTime <- 10  # last 10 years of data
  if        (tstep_pool == "weekly"){    timeFac <- 52
  } else if (tstep_pool == "timestep"){  timeFac <- 17520
  } else if (tstep_pool == "daily"){     timeFac <- 365
  } else if (tstep_pool == "monthly"){   timeFac <- 12
  } else if (tstep_pool == "yearly"){    timeFac <- 1 }
  spanRows <- timeFac * spanTime

  pdf("land_model_plots_prof.pdf", paper="a4", height=9.2, width=7.8, pointsize=12, 
      bg="white", family="Helvetica", version="1.6")
  print(paste0("Creating: land_model_plots_prof.pdf from the last ",spanTime," simulation years."))
  # loop over all output files
  for (i in 1:nOutFiles) {
    # only output files that contain layer information
    if(grepl("_lay",flist[i])) {
      print(paste('   reading',flist[i]))
      dat=read.table(flist[i],header=T,sep="")
      nvar=dim(dat)[2]
      # number of rows with data
      rowNum <- nrow(dat)
      # number of columns with variables, excluding the 1st row "Time"
      varCol <- ncol(dat) -1 
      # number of variables
      varNum <- varCol/ soilLayer
      # create one soil profile figure for each variable
      countVar <- 0
      for (v in 1:varNum) {
        a <- soilLayer *countVar +2
        b <- a + soilLayer -1
        vDat <- dat[,a:b]
        # calc start and end row to consider
        rowE <- nrow(vDat)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        # calc mean and SD per row for the last 10 years of simulated data, or all data im simulation time is < 10 yr
        vMean<-apply(vDat[rowS:rowE,], 2, mean)
        vSD  <-apply(vDat[rowS:rowE,], 2, sd)
        # prepare visualization of SD
        vSDup<-vMean+vSD
        vSDdn<-vMean-vSD
        # variable name
        vName<-colnames(vDat[1])
        vName<-substr(vName, 1, nchar(vName)-3)
        # set margins
        par(oma=c(0.2,0.2,0.2,0.2),mar=c(7.2,2.4,1.2,0.4),xpd=FALSE)
        # plot as line and points
        plot(vMean,soilLSeq,
             type="l",
             ylab="soil layers",xlab=vName,
             ylim=rev(range(soilLSeq)),
             xlim=c(min(vSDdn),max(vSDup)),
             main=paste0(vName,"    from ",flist[i]),
             tck=0.01, mgp=c(1.2, 0.14, 0))
        arrows(vSDup, soilLSeq, vSDdn, soilLSeq, code=3, length=0.0, angle=90, col="orange3")
        points(vMean,soilLSeq, pch=21)
        par(xpd=NA)
        legend(min(vSDdn), max(soilLSeq)+0.12*max(soilLSeq), 
                c(paste0("Mean over last ",spanTime," years."),"SD of mean",paste0("The mean was calculated from ",tstep_pool," model output.")), 
                cex=0.8, col=c("black","orange3",NA), pch=c(21,NA,NA), lty=c(NA,1,NA), box.col = "gray80")
        countVar <- countVar +1
      }
    }
  }
  garbage <- dev.off()
} # end: if (any(flistLayFiles))
