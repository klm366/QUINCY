## 
## calculate the partial correlation coefficient of the QUINCY parameters
##  from results of a parameter sensitivity simulation run
## 


## options ########################################################################################
# path to experiment directory (e.g. "/Net/Groups/BGI/scratch/USER/slm/exp_name/")
#expDir     <- "/Net/Groups/BGI/scratch/  /"   # must end with a slash ".../"  !
expDir     <- "/home/jengel/remote/rt-bgi-scratch-je/slm/exp_hainich_mavgparam_sensi_r1644/"   # must end with a slash ".../"  !
# output variable of which the partial correlation coefficient may be calculated (column name in output file)
outVarName <- "GPP"
# output file of output variabale; name without out-interval and ".txt" but must end with "_" examples: "vegfluxC_" or "soilfluxN_" 
outFile    <- "vegfluxC_"
# [T/F] run default set of output variables (overwriting the above definition of outVarName)
outVarSet  <- "F"
# years used for averaging the ouput value (at the end of the simulation results available)
avgYears   <- 10
# [T/F] transform both parameter values and output-variable values into ranks [highest value = highest rank] (bypassing potential problems from non-linearities) 
rankTform  <- "T"



## default set of output variables ################################################################
outVarNameSet <- c("GPP",
                   "Leaf_CN",
                   "Leaf_NP",
                   "Total_C",
                   "TotalSoilOrgC"
)

outFileSet    <- c("vegfluxC_",
                   "derivedLeafCN",   # the two input files needed here are hardcoded below
                   "derivedLeafNP",   # the two input files needed here are hardcoded below
                   "vegpoolC_",
                   "ecopoolCNP_"
)




                   
###################################################################################################
##          no need to change anything below this line                                           ##
###################################################################################################

# libraries
library(ppcor)

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# various global variables
outFileExt <- ".txt"



# change into the exp directory
setwd(expDir)


###################
# list directories (may cover all dirs from 1000_* to 1999_*)
dirList <- list.dirs(path = "./", full.names = F, recursive = F)
# remove "parameter_list" "x_binary" dirs from list
dirList <- dirList[1:(which(dirList == "parameter_list")-1)]
# get number of simulations
nSim    <- length(dirList)


###################
# get parameter_table.txt
paramTab    <- read.table("parameter_list/parameter_table.txt", header=T, sep="", stringsAsFactors=F)
# remove the first column from table with parameter values
paramTab    <- within(paramTab, rm("X0"))
# get parameter list
paramList   <- colnames(paramTab)


###################
# get output interval from namelist  (timestep daily weekly monthly yearly)
namelist  <- read.delim(paste0(dirList[1],"/input/namelist.slm"), header=F, sep="=", stringsAsFactors=F)
tstepPool <- trim(namelist[trim(namelist[,1])=="output_interval_pool",2])
# not differentiating pool vs flux tstep_* !


###################
# prepare reading the output-value data from all simulations into array

# create df for avg values of output variable
outVarAvgVal           <- data.frame(varName = numeric(), stringsAsFactors=F)
colnames(outVarAvgVal) <- outVarName[1]

# calc the number of lines from the output file that are used for averaging
 timeFac   <- 0    # scaling factor to get the correct number of rows
 if        (tstepPool == "weekly"){    timeFac <- 52
 } else if (tstepPool == "timestep"){  timeFac <- 17520
 } else if (tstepPool == "daily"){     timeFac <- 365
 } else if (tstepPool == "monthly"){   timeFac <- 12
 } else if (tstepPool == "yearly"){    timeFac <- 1 }
 spanRows  <- timeFac * avgYears



##############################
# loop over output variables #
##############################
nOutVar <- length(outVarName)
for (v in 1:nOutVar){
  

  ###################
  # create df for the data of all output variables


  ###################
  # store all the pcc results from the parameters in a list of df
  
  
  ###################
  # read the output-value data from all simulations into array and calc the avg value
  for (i in 1:nSim){
    print(paste0(" reading  ", dirList[i], "/output/", outFile, tstepPool, outFileExt, " <- ", outVarName[v]))
    dat <- read.table(paste0(dirList[i], "/output/", outFile, tstepPool, outFileExt), header=T, sep="")
    # get colnumber of outVarName[v]
    varColNumber <- which(colnames(dat) == outVarName[v])
    # test if enough data is available for defined avg-period
    rowE <- nrow(dat)
    if (rowE > spanRows){ rowS <- rowE - spanRows
    } else              { rowS <- 1  }
    # get avg value
    outVarAvgVal[i,1] <- mean(dat[rowS:rowE,varColNumber])
  }
  
   
  ###################
  # rank transform values of both parameters and output-variable
  outInfo <- "_no_ranktform"
  if (rankTform == "T"){
   print(" applying rank-transformation")
   outInfo <- "_with_ranktform"
   # transform output variable
   outVarAvgVal[,1] <- rank(outVarAvgVal[,1])
   # transform parameter values (sapply: apply function "rank" to each column/parameter)
   paramTab <- as.data.frame(sapply(paramTab, rank))
  }
    
  
  ###################
  # loop over parameters in paramList and calc the pcc for each of them
  
  # get numbers of parameters
  nParam <- length(paramList)
  
  # create results array
  pccResults           <- array(0, dim=c(nParam,3)) 
  colnames(pccResults) <- c("pcc_estimate", "p_value", "statistic")  # see ?pcor.test
  rownames(pccResults) <- colnames(paramTab)
  
  # run loop over parameters
  print(" running pcor.test() for each parameter")
  for (i in 1:nParam){
    # get only the values (as df) of the one parameter which the pcc is calculated for
    paramVal    <- subset(paramTab, select = paramList[i])
    # remove this parameter from table with parameter values
    paramTabRem <- paramTab[ , !(names(paramTab) %in% paramList[i])]
    # disable warnings of which pcor.test() creates a lot
    options(warn=-1)
    # run partial correlation analysis
    pcc <- pcor.test(outVarAvgVal[,1],
              paramVal[,1],
              paramTabRem[,1:ncol(paramTabRem)],
              method="pearson")
    # enable warnings
    options(warn=0)
    # save partial correlation coefficient in the results array
    pccResults[i,1] <- pcc[1,1]
    pccResults[i,2] <- pcc[1,2]
    pccResults[i,3] <- pcc[1,3]
    # note on pcc range:
    #  Correlation coefficients are always values between -1 and 1, 
    #  where -1 shows a perfect, linear negative correlation, and 1 shows a perfect, linear positive correlation.
  }

  # transform results to df
  pccResDF <- as.data.frame(pccResults)
  # sort with decreasing pcc_estimates
  pccResDF <- pccResDF[order(-pccResDF$pcc_estimate), ]
  

} ## end loop over output variables
###################################



###################
# visualisation and file output

# write pdf file
print(" create and write figure")
pdf(paste0("parameter_list/pcc_", outVarName, outInfo, ".pdf"), paper="a4r", height=6.4, width=12, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")
  par(mar=c(10.0,2.2,1.0,0.2) , oma=c(0.5,0.1,0.5,0.1)) # oma = outer margins (b, l, t, r)
  plot(1:nrow(pccResDF), pccResDF$pcc_estimate, 
     type="p", pch=15, cex=0.8,
     ylim=c(-1,1),
     main=paste0(outVarName, "   (* p<0.1)"),
     ylab=paste0("Part. Correl. Coeff:   X-parameter ~ avg ", outVarName),
     xlab=NA, #"QUINCY parameter",
     xaxt="n",
     mgp=c(1.2, 0.14, 0), tck=0.02
     )
  # plot a star for pcc_estimates < 0.1
  xAdd <- 0.1
  yAdd <- 0.04
  for (s in 1:nrow(pccResDF)){
    if (pccResDF$p_value[s] < 0.1){
      points(s+xAdd, pccResDF$pcc_estimate[s]+yAdd, pch=8, cex=0.6, col="dodgerblue3")
    }
  }
  #axis(1, at=1:nrow(pccResDF), labels=rownames(pccResDF))
  axis(1, at=seq(1, nrow(pccResDF), by=1), labels = FALSE, tck=0.02)
  text(seq(1, nrow(pccResDF), by=1), par("usr")[3] - 0.06, 
     labels=rownames(pccResDF), srt=65, pos=2, xpd = TRUE,
     cex=0.9)
  abline(h=c(-1.0,-0.5,0,0.5,1.0), col="sienna",lw=0.6)
garbage <- dev.off()

# write data file
print(" write output file")
pccResDFout              <- cbind(pccResDF,rownames(pccResDF))
colnames(pccResDFout)[4] <- "QUINCY_param"
write.table(pccResDFout, file=paste0("parameter_list/pcc_", outVarName, outInfo, ".csv"), append=F, 
            sep = ",", eol="\n", row.names=F, col.names=T)


print(" Done!")

# fin #
