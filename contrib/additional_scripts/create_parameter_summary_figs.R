##    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
##    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
##
##    This file is part of the QUINCY model.
##    The QUINCY model is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    The QUINCY model is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.

# 
# this script:
#     creates figures from your parameter sensitivity simulations
#     one file with a set of figures will be created per parameter in parameter_table.txt
# 
# run: from within the exp_xx/parameter_list/ folder 
#      Rscript --vanilla ./create_parameter_summary_figs.R
#
# note: Each figure shows the variation in results for a particular output variable covering 
#       the full parameter space created with the LHS and considering all parameters. 
#       The colors highlight the distribution of the values of the particular parameter across 
#       the range of output values of the variable, however, this does not necessarily indicate 
#       the effect of the parameter value on the variation of the output variable.
#


## define output file of which the variables should be used for fig creation
#  e.g.: vegpoolC_  or  vegfluxN_  or  soilpooldN15org_lay_  or  soilpoolPOrg_lay_  or  ecopoolCNP_  or ... I think you got the idea
outputFileName <- "vegpoolC_"

## show results of all simulations in the figure, set this variable to "T"
#  if set to "F" only the simulation with min and max values and about 8 values in between will be plotted
allSimulations <- "T"



#######################################################################################################################
## not change any below this line unless you know what your are doing                                                ##
#######################################################################################################################


## use for debugging and testing
#setwd("/home/jengel/remote/rt-bgi-scratch-je/slm/exp_testparamsensi_figscript/parameter_list")


## test if this script is executed from within the parameter_list/ dir with the file parameter_table.txt present
fileListTest <- Sys.glob("parameter_table.txt") 
if(!grepl("parameter_table.txt",fileListTest)) {
  print(" You need to run this script from within the folder parameter_list/ .")
  print(" EXIT!")
  exit("1")
}

## load library needed for color legend, using legend.gradient() function
library(SDMTools)

## create output folder for the figures
opath <- paste0(outputFileName,"sensi_figs")
dir.create(paste0("./",opath), showWarnings=T, recursive=T)

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

## read namelist to get output interval
namelistPath <- Sys.glob(paste0("../1001_*/input/namelist.slm"))
namelist     <- read.delim(namelistPath, header=F, sep="=", stringsAsFactors=F, skip=17)
tStepPool    <- trim(namelist[trim(namelist[,1])=="output_interval_pool",2])

## compile full output-file name
outputFileName <- paste0(outputFileName,tStepPool,".txt")
outputFilePath <- Sys.glob(paste0("../1001_*/output/",outputFileName))
outputFile     <- read.table(outputFilePath, header=T, sep="")
outputFileTime <- outputFile[,1]
nOutVar        <- ncol(outputFile) - 1
nOutValues     <- nrow(outputFile)

## message
print(paste0("Start sensitivity-simulation figure creation for variables from: ",outputFileName))
print(" ")

## read param table and extract information
paramTable              <- read.table("parameter_table.txt", header=T, sep="")
colnames(paramTable)[1] <- "simNo"
nParam                  <- ncol(paramTable) - 1
nSim                    <- nrow(paramTable)

## set number of simulations used for ploting lines of data
plotSim <- seq(1,nSim,1)
if (nSim > 12 & allSimulations != "T"){
  plotSim <- seq(1,nSim,2)
  plotSim <- c(plotSim[2:length(plotSim)],max(nSim))
}
if (nSim > 23 & allSimulations != "T"){
  plotSim <- seq(1,max(nSim),as.integer(nSim/12))
  plotSim <- c(plotSim[2:length(plotSim)],max(nSim))
}

## create colors for figs
colFunc <- colorRampPalette(c("red","yellow","blue2"))
colGrad <- colFunc(nSim)


## loop over parameters in the parameter_table.txt, create one figures-file per parameter
#  p is number of columns of the paramter
for (p in 2:(nParam+1)){
parameterName  <- colnames(paramTable)[p]
sortParamTable <- cbind(paramTable[,1],paramTable[,p]) 
sortParamTable <- sortParamTable[order(sortParamTable[,2]),] 

pdf(paste0(opath,"/",parameterName,".pdf"), paper="a4", height=8, width=9, pointsize=10, 
    bg="white", family="Helvetica", version="1.6")

## feedback on script progress
print(parameterName)

## loop over all output variables from the output file defined
#  v is number of columns of the variable
for (v in 2:(nOutVar+1)){
  
  ## create empty array for storing all values of the var from all simulations
  varData <- array(0, dim=c(nOutValues,nSim)) 
  
  ## get variables name
  varName  <- colnames(outputFile)[v]

  ## feedback on script progress
  print(paste0("  ",varName))

  ## loop over all simulations 1001 to 1xxx, but use sequence as in sortParamTable[,1] and gather values of the variables in an array
  for (s in 1:nSim){
    simNumber      <- 1000 + sortParamTable[s,1]
    outputFilePath <- Sys.glob(paste0("../",simNumber,"_*/output/",outputFileName))
    outputFile     <- read.table(outputFilePath, header=T, sep="")
    varData[,s]    <- outputFile[,v]
  }

  ## prepare figure creation
  yMinData <- min(varData)
  yMaxData <- max(varData)
  xMinData <- 0
  xMaxData <- nOutValues + (nOutValues * 0.09)  # add 9% to get space for legend
  xLegend1 <- nOutValues + (nOutValues * 0.025)
  xLegend2 <- nOutValues + (nOutValues * 0.055)
  yLegend1 <- yMinData + ((yMaxData-yMinData) * 0.3)
  yLegend2 <- yMaxData - ((yMaxData-yMinData) * 0.05)
  pnts     <- cbind(x = c(xLegend1,xLegend2,xLegend2,xLegend1), y =c(yLegend2,yLegend2,yLegend1,yLegend1))
    
  ## create figure
  plot(1:nOutValues, varData[,1], 
       type='l', lwd=0.6,
       col=colGrad[1],
       ylim=c(yMinData, yMaxData),
       xlim=c(xMinData, xMaxData),
       xlab=tStepPool, 
       ylab=colnames(outputFile)[v],
       main=paste0("Distribution of ", parameterName," values across the simulation results (see legend)"),
       tck=0.01, mgp=c(1.5, 0.2, 0))
  for (s in plotSim){
    lines(1:nOutValues, varData[,s], col=colGrad[s])
  }
  legend.gradient(pnts, colGrad, c(sortParamTable[1,2],sortParamTable[nSim,2]), title="values")
  
} # end loop output variables

trash <- dev.off()

} # end loop parameters


# finish message
print(" ")
print(" Done. Figures created.")
