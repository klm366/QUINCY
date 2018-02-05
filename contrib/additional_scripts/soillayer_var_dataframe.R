#!/usr/bin/env Rscript

# 
# this function does: create one big data.frame from all soil*[pool/flux]_lay_* files in location 'workDir'
# 



###############################################################################
# function start & change into 'workDir'                                      #
###############################################################################

getSoilLayVarDF <- function(workDir){

setwd(workDir)

  
###############################################################################
# get the output files and gather some information about the simulation run   #
###############################################################################

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# read timestep ans oil layer information of output file from namelist for pool & flux
namelist=read.delim('namelist.slm',header=F,sep="=",stringsAsFactors=F)
tstep_pool=trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
tstep_flux=trim(namelist[trim(namelist[,1])=="output_interval_flux",2])
soilLayer <- trim(namelist[trim(namelist[,1])=="nsoil_energy",2])
soilLayer <- as.numeric(soilLayer)
soilLSeq  <- seq(1,soilLayer,1)

# get soil*[pool/flux]_lay_* files
flist_pool = Sys.glob(paste0("soilpool*","_lay_",tstep_pool,"*.txt")) 
#flist_pool = c(flist_pool,Sys.glob(paste0("soilpool*","_",tstep_pool,"*.txt"))) 
# get flux output files  
flist_flux = Sys.glob(paste0("soilflux*","_lay_",tstep_flux,"*.txt"))
#flist_flux = c(flist_flux,Sys.glob(paste0("soilflux*","_",tstep_flux,"*.txt")))

# create overall file list
flist=c(flist_pool,flist_flux)
# number of output files
nOutFiles=length(flist)

if (nOutFiles < 1){
  print("no soil* lay files found") 
  EXIT 
}


###############################################################################
# read each individual file, transform and put into one data.frame            #
###############################################################################

print(" Gathered soil* lay file names. Start compiling the data.frame slVarDF.")
print("----------")

# create empty df
# columns: sl1 sl2 sl3 sl4 sl5 slx ... varID varName fileName time timeCont 
slVarDF <- data.frame()

# global variable ID starting with 1
varID <- 1

# loop over output files in the list
for (i in 1:nOutFiles) {
    print(paste(' working on ',flist[i]))
    dat=read.table(flist[i],header=T,sep="")
    # number of rows with data
    rowNum <- nrow(dat)
    # number of columns with variables, excluding the 1st row "Time"
    varCol <- ncol(dat) -1 
    # number of variables
    varNum <- varCol/ soilLayer
    # file name (removing all from _lay* on)
    fileName <- gsub("_lay.*", "", flist[i])

    # loop over variables in this file, only needed to get the column number right
    countVar <- 0
    for (v in 1:varNum) {
      # get column numbers for this variable
      a <- soilLayer *countVar +2
      b <- a + soilLayer -1
      # create extra df for the variable
      dat2 <- dat[,a:b]
      # variable name without soil layer information '_Lx' 
      vName<-colnames(dat2[1])
      vName<-substr(vName, 1, nchar(vName)-3)
      # add column names to variable df (simply numbers for sl numbers)
      colnames(dat2)<-soilLSeq

      # add columns to this dataframe of the variable
      # columns: varID varName fileName time timeCont
      dat2[,"varID"]    <- rep(varID,nrow(dat2))
      dat2[,"varName"]  <- rep(vName,nrow(dat2))
      dat2[,"filename"] <- rep(fileName,nrow(dat2))
      dat2[,"time"]     <- dat[,1]
      dat2[,"timeCont"] <- seq(dat[1,1],dat[1,1]*nrow(dat),dat[1,1])  # not very correct because of year-breaks
      
      # add this df to the global df
      slVarDF <- rbind(slVarDF,dat2)
      
      # increase variable-count      
      countVar <- countVar +1
      varID    <- varID    +1
    }
}

print("----------")
print(paste0(" DF slVarDF created, with ",varID," number of variables."))
print("----------")


return(slVarDF)
}

###############################################################################
# function end                                                                #
#                                                                             #
# start your own code here getting the slVarDF dataframe from the function    #
###############################################################################


# example function call:
#  exampleDF <- getSoilLayVarDF("/path/to/folder/with/results")

# you could also use this function in another R file like:
# 
#  source("/path/to/contrib/additional_scripts/soillayer_var_dataframe.R")
#  exampleDF <- getSoilLayVarDF("/path/to/folder/with/results")


