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

## 
## calculate the partial correlation coefficient of the QUINCY parameters
##  from results of a parameter sensitivity simulation run
## 

###################################################################################################
# unit conversion of output variables
molar_mass_C <- 12.0112
molar_mass_N <- 14.0072
molar_mass_P <- 30.997
weeks        <- 52
months       <- 12
outInterval  <- months    # update as needed
# soil layer depths (i.e. layer thickness) for 15 soil layers [m]
slDepths     <- c(0.065,0.065,0.065,0.069900094,0.09755344,0.136146793,0.190008156,0.265177744,0.370085354,0.516495718,0.720827842, 1.005996292,1.403980925,1.959413224,2.569414418)
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
                    "EcoC_1m",                  # 5  custom calc from vegPoolC and upper 1 m of soil C
                    "TotalSoilOrgC",            # 6
                    "Nfix_asymb",               # 7  asymbiontic N fixation only (SB var: SUM(asymb_n_fixation(ipoi,:)*soil_depth_sl(ipoi,:)))
                    "NetMinNH4",                # 8
                    "Net_Mineralisation_PO4",   # 9
                    "Nlosses",                  # 10  leaching + emmission | Emission_NOY Emission_N2O Emission_N2 LeachingNH4 LeachingNO3 LeachingDOM_N
                    "Plosses"                   # 11  leaching             | Leaching_PO4 LeachingDOM_P
)

outFileSet     <- c("vegfluxC_",                # 1
                    "derivedLeafCN",            # 2  vegpoolC_ & vegpoolN_  these two input files are hardcoded below
                    "derivedLeafNP",            # 3  vegpoolN_ & vegpoolP_  these two input files are hardcoded below
                    "vegpoolC_",                # 4
                    "ecopoolCNP_",              # 5
                    "ecopoolCNP_",              # 6
                    "soilfluxN_",               # 7
                    "soilfluxN_",               # 8
                    "soilfluxP_",               # 9
                    "soilfluxN_",               # 10  
                    "soilfluxP_"                # 11 
)

outVarUnitConvSet <- c(molar_mass_C / 1000 * outInterval,    # 1  mol/m2/week -> kg/m2/yr
                       molar_mass_C / molar_mass_N,          # 2  mol/mol -> g/g
                       molar_mass_N / molar_mass_P,          # 3  mol/mol -> g/g
                       molar_mass_C / 1000,                  # 4  mol C -> kg C
                       molar_mass_C / 1000,                  # 5  mol C -> kg C
                       molar_mass_C / 1000,                  # 6  mol C -> kg C
                       molar_mass_N * outInterval,           # 7  mol N/m2/week -> g/m2/yr
                       molar_mass_N * outInterval,           # 8  mol N/m2/week -> g/m2/yr
                       molar_mass_P * outInterval,           # 9  mol P/m2/week -> g/m2/yr
                       molar_mass_N * outInterval,           # 10 mol N/m2/week -> g/m2/yr
                       molar_mass_P * outInterval            # 11 mol P/m2/week -> g/m2/yr
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
namelist  <- read.delim(paste0(dirList[1],"/input/namelist.slm"), header=F, sep="=", stringsAsFactors=F, skip=17)
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
# create df for the pcc results only, and keep it unsortet (the pccResDFList keeps the sortet values)
pccResUsortDF   <- data.frame(varName = numeric(), stringsAsFactors=F)

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
  if (outVarName[v] == "Leaf_CN" || outVarName[v] == "Leaf_NP" || outVarName[v] == "EcoC_1m" || outVarName[v] == "Nlosses" || outVarName[v] == "Plosses"){
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
    # EcoC_1m
    if (outVarName[v] == "EcoC_1m"){
      outFileV      <- "vegpoolC_"
      outFileS      <- "soilpoolCOrg_lay_"
      outVarNameV   <- "Total_C"
      outVarNameSsl <- c("SoluLit_L1", "SoluLit_L2", "SoluLit_L3", "SoluLit_L4", "SoluLit_L5", "SoluLit_L6", "SoluLit_L7", "SoluLit_L8")
      outVarNameSpl <- c("PolyLit_L1", "PolyLit_L2", "PolyLit_L3", "PolyLit_L4", "PolyLit_L5", "PolyLit_L6", "PolyLit_L7", "PolyLit_L8")
      outVarNameSwl <- c("WodLit_L1",  "WodLit_L2",  "WodLit_L3",  "WodLit_L4",  "WodLit_L5",  "WodLit_L6",  "WodLit_L7",  "WodLit_L8")
      outVarNameSmi <- c("Mic_L1",     "Mic_L2",     "Mic_L3",     "Mic_L4",     "Mic_L5",     "Mic_L6",     "Mic_L7",     "Mic_L8")
      outVarNameSre <- c("Res_L1",     "Res_L2",     "Res_L3",     "Res_L4",     "Res_L5",     "Res_L6",     "Res_L7",     "Res_L8")
      outVarNameSmy <- c("Myc_L1",     "Myc_L2",     "Myc_L3",     "Myc_L4",     "Myc_L5",     "Myc_L6",     "Myc_L7",     "Myc_L8")
      soilCOutVarDF <- data.frame("solulit" = outVarNameSsl,
                                  "polylit" = outVarNameSpl,
                                  "wodlit"  = outVarNameSwl,
                                  "mic"     = outVarNameSmi,
                                  "res"     = outVarNameSre,
                                  "myc"     = outVarNameSmy,
                                  stringsAsFactors = F)
      print(paste0(" reading  ", dirList[1], "/output/", outFileV, tstepPool, outFileExt, " <- ", outVarNameV))
      print(paste0(" reading  ", dirList[1], "/output/", outFileS, tstepPool, outFileExt, " <- ", "sl1:8 (0.954 m) SoluLit PolyLit WodLit Mic Res Myc"))
      print(" ...")
      for (i in 1:nSim){
        # veg C
        dat1 <- read.table(paste0(dirList[i], "/output/", outFileV, tstepPool, outFileExt), header=T, sep="")
        varColNumber1 <- which(colnames(dat1) == outVarNameV)
        # soil C within upper 1 m (1:8 sl for 15 sl)
        dat2 <- read.table(paste0(dirList[i], "/output/", outFileS, tstepPool, outFileExt), header=T, sep="")
        # multiply columns by their depth (slDepths)
        nColDF        <- ncol(soilCOutVarDF)   # solulit polylit wodlit mic res myc
        nRowDF        <- nrow(soilCOutVarDF)   # L1 - L8
        varColNumber2 <- c(rep(0,nColDF*nRowDF))
        varCount      <- 0
        for (nc in 1:nColDF){
          for (nr in 1:nRowDF){
            varCount                       <- varCount + 1
            varColNumber2[varCount]        <- which(colnames(dat2) == soilCOutVarDF[nr,nc])#
            dat2[,varColNumber2[varCount]] <- dat2[,varColNumber2[varCount]] * slDepths[nr]
          }
        }
        rowE <- nrow(dat1)
        if (rowE > spanRows){ rowS <- rowE - spanRows
        } else              { rowS <- 1  }
        outVarAvgVal[i,1] <- mean(dat1[rowS:rowE,varColNumber1], na.rm=T)            # veg C
        for (m in varColNumber2){
          outVarAvgVal[i,1] <- outVarAvgVal[i,1] + mean(dat2[rowS:rowE,m], na.rm=T)  # add soil C
        }
      } # end  for()
    } # end EcoC_1m
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
  print(" running pcor.test() for each parameter  (with warnings disabled)")
  for (i in 1:nParam){
    # get only the values (as df) of the one parameter which the pcc is calculated for
    paramVal    <- subset(paramTab, select = paramList[i])
    # remove this parameter from table with parameter values
    paramTabRem <- paramTab[ , !(names(paramTab) %in% paramList[i])]
    # disable warnings of which pcor.test() creates a lot
    options(warn=-1)
    # run partial correlation analysis exclude all NA/NaN/Inf from the output values
    inclOutVal <- which(! is.na(outVarAvgVal[,1]) & ! is.nan(outVarAvgVal[,1]) & ! is.infinite(outVarAvgVal[,1]))
    pcc <- pcor.test(outVarAvgVal[inclOutVal,1],
              paramVal[inclOutVal,1],
              paramTabRem[inclOutVal,1:ncol(paramTabRem)],
              method=c("pearson", "kendall", "spearman"))
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

  # transform results to df and store results unsorted (needed for the summary pcc figure)
  pccResUsortDF[1:nParam,v]  <- as.data.frame(pccResults)[,1]
  colnames(pccResUsortDF)[v] <- outVarName[v]
  rownames(pccResUsortDF)    <- rownames(as.data.frame(pccResults))

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
# save out-variable mean values for later visualization, such as, histogram ....
outVarAvgValSetDF <- as.data.frame(outVarAvgValSet)



###################
# create and write figures
print(" create and write figures")
# output file-name identifier
outFileID <- outVarName[1]
if (outVarSet == "T"){
  outFileID <- "outVarSet"
}
# create pdf
pdf(paste0("parameter_list/pcc_", outFileID, outInfo, ".pdf"), paper="a4r", height=6.4, width=12, pointsize=8, 
    bg="white", family="Helvetica", version="1.6")
  par(mar=c(10.0,2.2,1.0,0.2) , oma=c(0.5,0.1,0.5,0.1)) # oma = outer margins (b, l, t, r)


  # plot boxplots  
  # mean of absolute output-variable values
  meanAbsVal <- sapply(outVarAvgValSet, mean, na.rm=T)
  # normalize the values using the mean, so that the mean of all boxplots is one and the variation is kept
  for (v in 1:nOutVar){
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
  abline(h=c(1.0), col="sienna", lw=0.6)  # mean value (is 1.0 because the values are normalized using the mean)
  abline(h=c(0.0,0.5,1.5,2.0), col="gray80", lw=0.4)  # mean value (is 1.0 because the values are normalized using the mean)
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


  # summary pcc estimate values
  # min-max span of normalized output-variable values
  spanAbsVal <- c(rep(NA,nOutVar))
  for (v in 1:nOutVar){
    spanAbsVal[v] <- abs(max(outVarAvgValSet[,v], na.rm=T) - min(outVarAvgValSet[,v], na.rm=T))
  }
  # weighted mean of each parameter oer all output variables (weighted by the min-max span of the the output-variable)
  pccWeightMean <- c(rep(NA,nParam))
  for (n in 1:nParam){
    pccWeightMean[n] <- weighted.mean(abs(pccResUsortDF[n,]), spanAbsVal)
  }
  # transform to df
  pccWeightMeanDF           <- as.data.frame(pccWeightMean)
  rownames(pccWeightMeanDF) <- rownames(pccResUsortDF)
  # and now sort (needs a second helper-column to work)
  pccWeightMeanDF              <- cbind.data.frame(pccWeightMeanDF,pccWeightMeanDF)
  colnames(pccWeightMeanDF)[2] <- "helper"
  pccWeightMeanDF              <- pccWeightMeanDF[order(-pccWeightMeanDF$pccWeightMean),]
  # plot
  plot(1:nrow(pccWeightMeanDF), pccWeightMeanDF$pccWeightMean,  
       type="p", pch=15, cex=0.8,
       ylim=c(0,1),
       main="weighted mean of absolute! parameter pcc estimates (across output variables)",
       ylab="Part. Correl. Coeff of parameters (weighted mean)",
       xlab=NA, #"QUINCY parameter",
       xaxt="n",
       mgp=c(1.2, 0.14, 0), tck=0.02
  )
  axis(1, at=seq(1, nrow(pccWeightMeanDF), by=1), labels = FALSE, tck=0.02)
  text(seq(1, nrow(pccWeightMeanDF), by=1), par("usr")[3] - 0.06, 
       labels=rownames(pccWeightMeanDF), srt=65, pos=2, xpd = TRUE,
       cex=1.0)
  abline(h=c(0,0.25,0.5,0.75,1.0), col="sienna",lw=0.6)

  
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
      if (! is.nan(pccResDF$p_value[s]) & ! is.infinite(pccResDF$p_value[s])){
        if (pccResDF$p_value[s] < 0.1){
          points(s+xAdd, pccResDF$pcc_estimate[s]+yAdd, pch=8, cex=0.6, col="dodgerblue3")
        }
      }
    }
    #axis(1, at=1:nrow(pccResDF), labels=rownames(pccResDF))
    axis(1, at=seq(1, nrow(pccResDF), by=1), labels = FALSE, tck=0.02)
    text(seq(1, nrow(pccResDF), by=1), par("usr")[3] - 0.06, 
       labels=rownames(pccResDF), srt=65, pos=2, xpd = TRUE,
       cex=1.0)
    abline(h=c(-1.0,-0.5,0,0.5,1.0), col="sienna",lw=0.6)
  } # end loop over output variables


  # loop over output variables to plot histogram of the mean-values from all simulations
  for (v in 1:nOutVar){
    minOutVal    <- min(outVarAvgValSetDF[,v], na.rm=T)
    maxOutVal    <- max(outVarAvgValSetDF[,v], na.rm=T)
    meanOutVal   <- mean(outVarAvgValSetDF[,v], na.rm=T)
    medOutVal    <- median(outVarAvgValSetDF[,v], na.rm=T)
    oVarNameHist <- colnames(outVarAvgValSetDF)[v]
    hist(outVarAvgValSetDF[,v], 
         breaks=100,
         xlim=c(minOutVal,maxOutVal),
         main=oVarNameHist, xlab="mean of absolute values from individual parameter-set simulations")
    abline(v=c(minOutVal,maxOutVal,meanOutVal), col="sienna",lw=0.6)
    abline(v=medOutVal, col="dodgerblue3",lw=1.0)
    legend("topleft", inset=0.01, c("min, max, mean","median"), lty=c(1,1), col=c("sienna","dodgerblue3"), box.col = "gray84")
  }
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
# write normalized values of output variables (normalized by the mean)
write.table(outVarAvgValSet, file=paste0("parameter_list/out_variables_values_normalized_by_mean", outInfo, ".csv"),
            append=F, sep = ",", eol="\n", row.names=F, col.names=T)
# write absolute values of output variables per simulation (i.e., mean over last x years per output-value from each simulation)
write.table(outVarAvgValSetDF, file=paste0("parameter_list/out_variables_values_absolute", outInfo, ".csv"),
            append=F, sep = ",", eol="\n", row.names=F, col.names=T)            
# write absolute mean values of output variables
rowNames <- colnames(outVarAvgValSet)
write.table(meanAbsVal, file=paste0("parameter_list/out_variables_absolute_values_mean", outInfo, ".csv"),
            append=F, sep = ",", eol="\n", row.names=rowNames, col.names=F)
# write the values of "weighted mean of absolute! parameter pcc estimates (across output variables)" (pccWeightMeanDF)
write.table(pccWeightMeanDF$pccWeightMean, file=paste0("parameter_list/pcc_weighted_mean_all_parameters", outInfo, ".csv"), 
            append=F, sep = ",", eol="\n",
            row.names=rownames(pccWeightMeanDF), 
            col.names=F)



print(" Done!")

# fin #
