## 
## calculate the partial correlation coefficient of the QUINCY parameters
##  from results of a parameter sensitivity simulation run
## 

###################################################################################################
# helpers for unit conversion of output variables
molar_mass_C <- 12.0112
molar_mass_N <- 14.0072
molar_mass_P <- 30.997
weeks        <- 52
months       <- 12
outInterval  <- weeks    # update as needed
###################################################################################################



## options ########################################################################################
# path to experiment directory (e.g. "/Net/Groups/BGI/scratch/USER/slm/exp_name/")
expDir     <- "/Net/Groups/BGI/scratch/USER...  /"   # must end with a slash ".../"  !
# output variable of which the partial correlation coefficient may be calculated (column name in output file)
outVarName     <- "GPP"
# output file of output variabale; name without out-interval and ".txt" but must end with "_" examples: "vegfluxC_" or "soilfluxN_" 
outFile        <- "vegfluxC_"
# factor for unit conversion of the output variable
outVarUnitConv <- molar_mass_C / 1000 * outInterval
# [T/F] run default set of output variables (overwriting the above definition of outVarName)
outVarSet      <- "T"
# years used for averaging the ouput value (at the end of the simulation results available)
avgYears       <- 10
# [T/F] transform both parameter values and output-variable values into ranks [highest value = highest rank] (bypassing potential problems from non-linearities) 
rankTform      <- "T"



## default set of output variables ################################################################
outVarNameSet  <- c("GPP",                      # 1
                    "Leaf_CN",                  # 2  Leaf_C / Leaf_N
                    "Leaf_NP",                  # 3  Leaf_N / Leaf_P
                    "Total_C",                  # 4
                    "TotalSoilOrgC",            # 5
#                    "FixationN",                # 6  symbiontic N fixation (variable veg_mem n_fixation)
                    "NetMinNH4",                # 7
                    "Net_Mineralisation_PO4",   # 8
                    "Nlosses",                  # 9   leaching + emmission | Emission_NOY Emission_N2O Emission_N2 LeachingNH4 LeachingNO3 LeachingDOM_N
                    "Plosses"                   # 10  leaching             | Leaching_PO4 LeachingDOM_P
)

outFileSet     <- c("vegfluxC_",                # 1
                    "derivedLeafCN",            # 2  vegpoolC_ & vegpoolN_  these two input files are hardcoded below
                    "derivedLeafNP",            # 3  vegpoolN_ & vegpoolP_  these two input files are hardcoded below
                    "vegpoolC_",                # 4
                    "ecopoolCNP_",              # 5
#                    "vegfluxN_",                # 6
                    "soilfluxN_",               # 7
                    "soilfluxP_",               # 8
                    "soilfluxN_",               # 9  
                    "soilfluxP_"                # 10 
)

outVarUnitConvSet <- c(molar_mass_C / 1000 * outInterval,    # 1  mol/m2/week -> kg/m2/yr
                       molar_mass_C / molar_mass_N,          # 2
                       molar_mass_N / molar_mass_P,          # 3
                       molar_mass_C / 1000,                  # 4  mol C -> kg C
                       molar_mass_C / 1000,                  # 5  mol C -> kg C
#                       molar_mass_N * outInterval,           # 6  mol N/m2/week -> g/m2/yr
                       molar_mass_N * outInterval,           # 7  mol N/m2/week -> g/m2/yr
                       molar_mass_P * outInterval,           # 8  mol P/m2/week -> g/m2/yr
                       molar_mass_N * outInterval,           # 9  mol N/m2/week -> g/m2/yr
                       molar_mass_P * outInterval            # 10 mol P/m2/week -> g/m2/yr
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

# use output-variable set if outVarSet=T
if (outVarSet == "T"){
  outVarName     <- outVarNameSet
  outFile        <- outFileSet
  outVarUnitConv <- outVarUnitConvSet
}
  
# number of output variables to cover 
nOutVar <- length(outVarName)

###################
# create df for the data of all output variables
outVarAvgValSet <- data.frame(varName = numeric(), stringsAsFactors=F)

###################
# store all the pcc results from the parameters in a list of df
pccResDFList <- list()

###################
# rank transform values of parameters already
outInfo <- "_no_ranktform"
if (rankTform == "T"){
  print(" applying rank-transformation to parameter values")
  outInfo <- "_with_ranktform"
  # transform parameter values (sapply: apply function "rank" to each column/parameter)
  paramTab <- as.data.frame(sapply(paramTab, rank))
}


# loop over output variables
print(" reading the output-variable data from all simulations into array and calc the avg value of each")
for (v in 1:nOutVar){

  # create df for avg values of output variable
  outVarAvgVal           <- data.frame(varName = numeric(), stringsAsFactors=F)
  colnames(outVarAvgVal) <- outVarName[v]
  
  ###################
  # read the output-variable data from all simulations into array and calc the avg value of each
  # diff normal & derived variables (case 1 and 2)
  derivedOutVar <- "F"
  if (outVarName[v] == "Leaf_CN" || outVarName[v] == "Leaf_NP" || outVarName[v] == "Nlosses" || outVarName[v] == "Plosses"){
    derivedOutVar <- "T"
  }
  
  # case 1: one variable (default case)
  if (derivedOutVar == "F"){
    print(paste0(" reading  ", dirList[1], "/output/", outFile[v], tstepPool, outFileExt, " <- ", outVarName[v]))
    print(" ...")
    for (i in 1:nSim){
      #print(paste0(" reading  ", dirList[i], "/output/", outFile[v], tstepPool, outFileExt, " <- ", outVarName[v]))
      dat <- read.table(paste0(dirList[i], "/output/", outFile[v], tstepPool, outFileExt), header=T, sep="")
      # get colnumber of outVarName[v]
      varColNumber <- which(colnames(dat) == outVarName[v])
      # test if enough data is available for defined avg-period
      rowE <- nrow(dat)
      if (rowE > spanRows){ rowS <- rowE - spanRows
      } else              { rowS <- 1  }
      # get avg value
      outVarAvgVal[i,1] <- mean(dat[rowS:rowE,varColNumber], na.rm=T)
    } # end  for()
  } # end  if()

  # case 2: derived variable from more than one variables (and output files) ####################################################
  if (derivedOutVar == "T"){                                                                                                   ##
    # Leaf_CN
    if (outVarName[v] == "Leaf_CN"){
      outFileC    <- "vegpoolC_"
      outFileN    <- "vegpoolN_"
      outVarNameC <- "Leaf_C"
      outVarNameN <- "Leaf_N"
      print(paste0(" reading  ", dirList[1], "/output/", outFileC, tstepPool, outFileExt, " <- ", outVarNameC))
      print(paste0(" reading  ", dirList[1], "/output/", outFileN, tstepPool, outFileExt, " <- ", outVarNameN))
      print(" ...")
      for (i in 1:nSim){
        dat1 <- read.table(paste0(dirList[i], "/output/", outFileC, tstepPool, outFileExt), header=T, sep="")
        varColNumber1 <- which(colnames(dat1) == outVarNameC)
        dat2 <- read.table(paste0(dirList[i], "/output/", outFileN, tstepPool, outFileExt), header=T, sep="")
        varColNumber2 <- which(colnames(dat2) == outVarNameN)
        rowE <- nrow(dat1)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        outVarAvgVal[i,1] <- mean(dat1[rowS:rowE,varColNumber1], na.rm=T) / mean(dat2[rowS:rowE,varColNumber2], na.rm=T)  # C/N
      } # end  for()
    } # end Leaf_CN
    # Leaf_NP
    if (outVarName[v] == "Leaf_NP"){
      outFileN    <- "vegpoolN_"
      outFileP    <- "vegpoolP_"
      outVarNameN <- "Leaf_N"
      outVarNameP <- "Leaf_P"
      print(paste0(" reading  ", dirList[1], "/output/", outFileN, tstepPool, outFileExt, " <- ", outVarNameN))
      print(paste0(" reading  ", dirList[1], "/output/", outFileP, tstepPool, outFileExt, " <- ", outVarNameP))
      print(" ...")
      for (i in 1:nSim){
        dat1 <- read.table(paste0(dirList[i], "/output/", outFileN, tstepPool, outFileExt), header=T, sep="")
        varColNumber1 <- which(colnames(dat1) == outVarNameN)
        dat2 <- read.table(paste0(dirList[i], "/output/", outFileP, tstepPool, outFileExt), header=T, sep="")
        varColNumber2 <- which(colnames(dat2) == outVarNameP)
        rowE <- nrow(dat1)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        outVarAvgVal[i,1] <- mean(dat1[rowS:rowE,varColNumber1], na.rm=T) / mean(dat2[rowS:rowE,varColNumber2], na.rm=T)  # N/P
      } # end  for()
    } # end Leaf_NP
    # Nlosses
    if (outVarName[v] == "Nlosses"){
      outVarNames  <- c("Emission_NOY", "Emission_N2O", "Emission_N2", "LeachingNH4", "LeachingNO3", "LeachingDOM_N")
      varColNumber <- c(rep(0,length(outVarNames)))
      print(paste0(" reading  ", dirList[1], "/output/", outFile[v], tstepPool, outFileExt, " <- ", outVarNames))
      print(" ...")
      for (i in 1:nSim){
        dat <- read.table(paste0(dirList[i], "/output/", outFile[v], tstepPool, outFileExt), header=T, sep="")
        varColNumber[1] <- which(colnames(dat) == outVarNames[1])
        varColNumber[2] <- which(colnames(dat) == outVarNames[2])
        varColNumber[3] <- which(colnames(dat) == outVarNames[3])
        varColNumber[4] <- which(colnames(dat) == outVarNames[4])
        varColNumber[5] <- which(colnames(dat) == outVarNames[5])
        varColNumber[6] <- which(colnames(dat) == outVarNames[6])
        rowE <- nrow(dat)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        outVarAvgVal[i,1] <- mean(dat[rowS:rowE,varColNumber[1]], na.rm=T) + 
                             mean(dat[rowS:rowE,varColNumber[2]], na.rm=T) +
                             mean(dat[rowS:rowE,varColNumber[3]], na.rm=T) +
                             mean(dat[rowS:rowE,varColNumber[4]], na.rm=T) +
                             mean(dat[rowS:rowE,varColNumber[5]], na.rm=T) +
                             mean(dat[rowS:rowE,varColNumber[6]], na.rm=T)
      } # end  for()
    } # end Nlosses
    # Plosses
    if (outVarName[v] == "Plosses"){
      outVarNames  <- c("Leaching_PO4", "LeachingDOM_P")
      varColNumber <- c(rep(0,length(outVarNames)))
      print(paste0(" reading  ", dirList[1], "/output/", outFile[v], tstepPool, outFileExt, " <- ", outVarNames))
      print(" ...")
      for (i in 1:nSim){
        dat <- read.table(paste0(dirList[i], "/output/", outFile[v], tstepPool, outFileExt), header=T, sep="")
        varColNumber[1] <- which(colnames(dat) == outVarNames[1])
        varColNumber[2] <- which(colnames(dat) == outVarNames[2])
        rowE <- nrow(dat)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        outVarAvgVal[i,1] <- mean(dat[rowS:rowE,varColNumber[1]], na.rm=T) + 
                             mean(dat[rowS:rowE,varColNumber[2]], na.rm=T)
      } # end  for()
    } # end Plosses
  } # end  'if (derivedOutVar == "T"){'                                                                                        ##
  # end case 2: derived variable from more than one variables (and output files) ################################################
  
  # add output-variable avg values to df befor rank-transformation is done
  outVarAvgValSet[1:nSim,v]    <- outVarAvgVal[,1]
  colnames(outVarAvgValSet)[v] <- colnames(outVarAvgVal)[1]
    
   
  ###################
  # rank transform also the values of avg output-variable
  if (rankTform == "T"){
   print(" applying rank-transformation to avg values of output variables")
   # transform output variable
   outVarAvgVal[,1] <- rank(outVarAvgVal[,1])
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

  # add this df to the list of df
  pccResDFList[[v]] <- pccResDF
    
  print("                                       ")
} ## end loop over output variables
###################################




###################
# unit conversion of output-variable values  (generally mol to gramm by multiplication with molar_massX)
print(" applying unit conversion to output variables")
for (v in 1:nOutVar){
  outVarAvgValSet[,v] <- outVarAvgValSet[,v] * outVarUnitConv[v]
}


###################
# create and write figures
print(" create and write figures")
# output file-name identifier
outFileID <- outVarName[1]
if (outVarSet == "T"){
  outFileID <- "outVarSet"
}
# create pdf
pdf(paste0("parameter_list/pcc_", outFileID, outInfo, ".pdf"), paper="a4r", height=6.4, width=12, pointsize=10, 
    bg="white", family="Helvetica", version="1.6")
  par(mar=c(10.0,2.2,1.0,0.2) , oma=c(0.5,0.1,0.5,0.1)) # oma = outer margins (b, l, t, r)

  # plot boxplots
  meanAbsVal <- sapply(outVarAvgValSet, mean, na.rm=T)
  for (v in 1:nOutVar){
    # normalize the values using the mean, so that the mean of all boxplots is one and the variation is kept
    outVarAvgValSet[,v] <- outVarAvgValSet[,v] / mean(outVarAvgValSet[,v],na.rm=T)
  }
  #minRelVal <- min(outVarAvgValSet)
  #maxRelVal <- max(outVarAvgValSet)
  boxplot(outVarAvgValSet, 
          #ylim=c(minRelVal-abs(0.0006*minRelVal),maxRelVal+abs(0.0006*maxRelVal)),
          main="Variation in output variables (each normalized by the mean)",
          col="gray80",
          range=0,                      # whiskers cover min and max value
          xlab=" \n Boxplots show: median, 25th and 75th percentile, min/max values \n mean = line at 1.0, and absolute mean values at x-axis",
          xaxt="n")
  abline(h=1.0, col="sienna", lw=0.6)  # mean value (is 1.0 because teh values are normalized using the mean)
  # x-axis labels: output-variable names & absolute mean values
  axis(1,at=c(1:ncol(outVarAvgValSet)), labels=colnames(outVarAvgValSet), tick=F, tck=0.02, mgp=c(2,0.2,0))
  # disable warnings to suppress the warning about mpg() using negative and positve values
  options(warn=-1)
  axis(1,at=c(1:ncol(outVarAvgValSet)), labels=round(meanAbsVal,2), tick=F, tck=0.02, mgp=c(2,-1.0,0))
  # enable warnings
  options(warn=0)
  # another way to plot the absolute mean values
  #text(seq(1, nrow(pccResDF), by=1), 
  #     par("usr")[3] + abs(0.0005*minRelVal), 
  #     labels=round(meanAbsVal,2),
  #     srt=0, pos=1, xpd = TRUE, cex=0.9)


  # loop over output variables to plot pcc_estimate values
  for (v in 1:nOutVar){
    # get df
    pccResDF <- pccResDFList[[v]]
    # plot pcc_estimate values
    plot(1:nrow(pccResDF), pccResDF$pcc_estimate, 
       type="p", pch=15, cex=0.8,
       ylim=c(-1,1),
       main=paste0(outVarName[v], "   (* p<0.1)"),
       ylab=paste0("Part. Correl. Coeff:   X-parameter ~ avg ", outVarName[v]),
       xlab=NA, #"QUINCY parameter",
       xaxt="n",
       mgp=c(1.2, 0.14, 0), tck=0.02
       )
    # plot a star for pcc_estimates < 0.1
    xAdd <- 0.1
    yAdd <- 0.04
    for (s in 1:nrow(pccResDF)){
      if (! is.nan(pccResDF$p_value[s]) && ! is.infinite(pccResDF$p_value[s])){
        if (pccResDF$p_value[s] < 0.1){
          points(s+xAdd, pccResDF$pcc_estimate[s]+yAdd, pch=8, cex=0.6, col="dodgerblue3")
        }
      }
    }
    #axis(1, at=1:nrow(pccResDF), labels=rownames(pccResDF))
    axis(1, at=seq(1, nrow(pccResDF), by=1), labels = FALSE, tck=0.02)
    text(seq(1, nrow(pccResDF), by=1), par("usr")[3] - 0.06, 
       labels=rownames(pccResDF), srt=65, pos=2, xpd = TRUE,
       cex=0.9)
    abline(h=c(-1.0,-0.5,0,0.5,1.0), col="sienna",lw=0.6)
  } # end loop over output variables
garbage <- dev.off()


###################
# data file output
print(" write output files")
# loop over output variables
for (v in 1:nOutVar){
  # get df
  pccResDF <- pccResDFList[[v]]
  # add new column with parameer names
  pccResDFout              <- cbind(pccResDF,rownames(pccResDF))
  colnames(pccResDFout)[4] <- "QUINCY_param"
  write.table(pccResDFout, file=paste0("parameter_list/pcc_", outVarName[v], outInfo, ".csv"), append=F, 
            sep = ",", eol="\n", row.names=F, col.names=T)
}





print(" Done!")

# fin #
