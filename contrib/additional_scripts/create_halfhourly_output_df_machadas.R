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

# 
# this script creates an R array with 1/2 hourly data from default QUINCY land model output with period timestep
#   input:  default QUINCY land model output with period timestep
#           range of years does not matter
#   output: array written into a .csv file as table with 1/2 hourly model output
#           filename: quincy_hh_timestep_land_output_[date]_for_machadas.csv
#           output variables are selected for comparison with field data measures at the Majadas site
#           units of model output are converted to match the units of field data
#           also creating a .pdf file with a figure for each variable (quincy_hh_timestep_land_output_[date]_for_machadas.csv.pdf)
#
# run this script:
#     copy this script into the directory with the QUINCY model output (i.e., containing files like *_timestep.txt)
#     cd into this directory
#     Rscript --vanilla ./create_halfhourly_output_df_machadas
#



## main variables

# input file extension (i.e., the quincy default output files used as input by this script)
inFileExt <- "_timestep.txt"
# time stamp
date <- Sys.Date()
# output file
outFile <- paste0("./quincy_hh_timestep_land_output_", date, "_for_machadas.csv")
# list with units of all output variables (is created below)
unitList <- "empty"


## soil layer width [m] (15 sl, 9.5 m depth)
soilLayWidth <- c(0.065, 0.065, 0.065, 0.069900094, 0.09755344, 0.136146793, 0.190008156, 0.265177744, 0.370085354, 0.516495718, 
                  0.720827842, 1.005996292, 1.403980925, 1.959413224, 2.569414418) 

## constants for unit conversion
molar_mass_C <- 12.0112
molar_mass_N <- 14.0072
molar_mass_P <- 30.997


## create output array

# number of soil layers
sPoolCOrg <- read.table(paste0("soilpoolCOrg_lay", inFileExt), header=T, sep="")
nSoilLay  <- (length(colnames(sPoolCOrg)) -1) / 10  # time + 10 variables: SoluLit PolyLit WodLit DOM aDOM Fungi Myc Mic Res aRes
print(paste0("number of soil layers is: ", nSoilLay))

# number of columns of the output array (3D variables store the layer-data in multiple columns)
nOutCol  <- 250   # this includes overprovisioning, columns without data will be removed before write the output file (see line ~370)
print(paste0("using ", nOutCol, " columns for the output array, including overprovisioning"))
# number of rows of the output array (i.e., number of timesteps from output data)
nOutRow  <- nrow(sPoolCOrg)
remove(sPoolCOrg)
# array
outArr   <- array(NA, dim=c(nOutRow, nOutCol))
# column names with just numbers for now (replaced by variable names below)
colnames(outArr) <- seq(1,nOutCol,1)
# column-count variable (staring with 1; the 1st var goes into the 1st column)
countCol <- 1


## get canopy layers (for calc of chlorophyll)
vegCanopy   <- read.table(paste0("veg_canopy", inFileExt), header=T, sep="")
nCanopyLay  <- nrow(vegCanopy) / nOutRow                                        # should be 10 by default
print(paste0("number of canopy layers is: ", nCanopyLay))
remove(vegCanopy)


## define function to get data of a specific variable
getVarData <- function(quincyOutputFile, variable){
  print(paste0(" reading ", variable, " from quincy output file ", quincyOutputFile))
  dat     <- read.table(quincyOutputFile, header=T, sep="")
  nColumn <- which(colnames(dat) == variable)
  if (length(nColumn) != 1){
    if (length(nColumn) < 1) stop("no such variable present in the file, check header names")
    if (length(nColumn) > 1) stop("more than two variables with the same name are present in the file, check header names")
  }
  varData <- dat[,nColumn]
  return(varData)
}


###
### go through the QUINCY default-output files and put variables into the outArr
###

##
inFile  <- paste0("soilfluxCdC13DC14", inFileExt)
var1    <- getVarData(inFile, "HetResp")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "Rh"               # µmol m-2 s-1
unitList                    <- "µmol m-2 s-1"     # initial entry of unitList object
countCol                    <- countCol + 1

##
inFile  <- paste0("soilpoolCOrg_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("SoluLit_L", l))
  var2    <- getVarData(inFile, paste0("PolyLit_L", l))
  var3    <- getVarData(inFile, paste0("WodLit_L", l))
  print(paste0("calc and save Litter_C for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2 + var3
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-3 -> kg C m-3
  colnames(outArr)[countCol]  <- paste0("Litter_C_lay", l)
  unitList                    <- c(unitList, "kg C m-3")
  countCol                    <- countCol + 1
}
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("Mic_L", l))
  var2    <- getVarData(inFile, paste0("Res_L", l))
  print(paste0("calc and save SOM_C for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-3 -> kg C m-3
  colnames(outArr)[countCol]  <- paste0("SOM_C_lay", l)
  unitList                    <- c(unitList, "kg C m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolNOrg_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("SoluLit_L", l))
  var2    <- getVarData(inFile, paste0("PolyLit_L", l))
  var3    <- getVarData(inFile, paste0("WodLit_L", l))
  print(paste0("calc and save Litter_N for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2 + var3
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_N / 1000   # mol N m-3 -> kg N m-3
  colnames(outArr)[countCol]  <- paste0("Litter_N_lay", l)
  unitList                    <- c(unitList, "kg N m-3")
  countCol                    <- countCol + 1
}
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("Mic_L", l))
  var2    <- getVarData(inFile, paste0("Res_L", l))
  print(paste0("calc and save SOM_N for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_N / 1000   # mol N m-3 -> kg N m-3
  colnames(outArr)[countCol]  <- paste0("SOM_N_lay", l)
  unitList                    <- c(unitList, "kg N m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolPOrg_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("SoluLit_L", l))
  var2    <- getVarData(inFile, paste0("PolyLit_L", l))
  var3    <- getVarData(inFile, paste0("WodLit_L", l))
  print(paste0("calc and save Litter_P for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2 + var3
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_P / 1000   # mol P m-3 -> kg P m-3
  colnames(outArr)[countCol]  <- paste0("Litter_P_lay", l)
  unitList                    <- c(unitList, "kg P m-3")
  countCol                    <- countCol + 1
}
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("Mic_L", l))
  var2    <- getVarData(inFile, paste0("Res_L", l))
  print(paste0("calc and save SOM_P for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_P / 1000   # mol P m-3 -> kg P m-3
  colnames(outArr)[countCol]  <- paste0("SOM_P_lay", l)
  unitList                    <- c(unitList, "kg P m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolInNOrg_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("NH4solute_L", l))
  var2    <- getVarData(inFile, paste0("NO3solute_L", l))
  print(paste0("calc and save Nmineral for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_N / 1000   # mol N m-3 -> kg N m-3
  colnames(outArr)[countCol]  <- paste0("Nmineral_lay", l)
  unitList                    <- c(unitList, "kg N m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolInPOrg_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("PO4assoc_fast_L", l))
  var2    <- getVarData(inFile, paste0("PO4assoc_slow_L", l))
  var3    <- getVarData(inFile, paste0("PO4solute_L", l))
  print(paste0("calc and save OlsonP for soil layer: ", l))
  outArr[,countCol]           <- var1 + var2 + var3
  outArr[,countCol]           <- outArr[,countCol] * molar_mass_P / 1000   # mol P m-3 -> kg P m-3
  colnames(outArr)[countCol]  <- paste0("OlsonP_lay", l)
  unitList                    <- c(unitList, "kg P m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolPhysics_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("WaterContentSoil_L", l))
  outArr[,countCol]           <- var1 / soilLayWidth[l]       # convert from meter to volumetric fraction
  colnames(outArr)[countCol]  <- paste0("Theta_lay", l)       # volumetric fraction
  unitList                    <- c(unitList, "vol. fraction")
  countCol                    <- countCol + 1
}
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("TemperatureSoil_L", l))
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- paste0("SoilT_lay", l)       # K
  unitList                    <- c(unitList, "K")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("soilpoolPhysics", inFileExt)
var1    <- getVarData(inFile, "RootZoneSoilWaterPotential")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "RootZone_SoilWater_Potential" # MPa
unitList                    <- c(unitList, "MPa")
countCol                    <- countCol + 1

##
inFile  <- paste0("soil_diagnostics_lay", inFileExt)
for (l in 1:nSoilLay){
  var1    <- getVarData(inFile, paste0("RootC_L", l))         # quincy output: veg_pool%fine_root%carbon(ipoi) * root_fraction_sl(ipoi,:) / soil_depth_sl(ipoi,:)
  outArr[,countCol]           <- var1 * molar_mass_C / 1000   # mol C m-3 -> kg C m-3
  colnames(outArr)[countCol]  <- paste0("RootDensity_lay", l)
  unitList                    <- c(unitList, "kg C m-3")
  countCol                    <- countCol + 1
}

##
inFile  <- paste0("veg_diagnostics", inFileExt)
var1    <- getVarData(inFile, "Height")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "Height"           # m
unitList                    <- c(unitList, "m")
countCol                    <- countCol + 1

var1    <- getVarData(inFile, "Diameter")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "Diameter"         # m
unitList                    <- c(unitList, "m")
countCol                    <- countCol + 1

var1    <- getVarData(inFile, "LAI")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "LAI"              # m2 m-2
unitList                    <- c(unitList, "m2 m-2")
countCol                    <- countCol + 1

##
inFile  <- paste0("vegfluxC", inFileExt)
var1    <- getVarData(inFile, "GPP")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "GPP"              # µmol m-2 s-1
unitList                    <- c(unitList, "µmol m-2 s-1")
countCol                    <- countCol + 1

var1    <- getVarData(inFile, "MaintResp")
var2    <- getVarData(inFile, "NTransformResp")
var3    <- getVarData(inFile, "GrowthResp")
outArr[,countCol]           <- var1 + var2 + var3
colnames(outArr)[countCol]  <- "Ra"               # µmol m-2 s-1
unitList                    <- c(unitList, "µmol m-2 s-1")
countCol                    <- countCol + 1

var1    <- getVarData(inFile, "Growth_sapwood_C")
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "NPP_sapwood"      # µmol m-2 s-1
unitList                    <- c(unitList, "µmol m-2 s-1")
countCol                    <- countCol + 1

##
inFile  <- paste0("vegfluxEnergy", inFileExt)
varList     <- c("Rnet", "Qh", "Qle", "Qg", "Albedo", "PPFD", "aPPFD", "fAPPFD", "albedo_vis", "albedo_nir")
varOutList  <- c("Rnet", "Qh", "Qle", "Qg", "Albedo", "downward_PPFD", "absorbed_PPFD", "fPAR", "albedo_vis", "albedo_nir")
unitVegFlux <- c("Wm-2", "Wm-2", "Wm-2", "Wm-2", "--", "µmol m-2 s-1", "µmol m-2 s-1", "--", "--", "--")
for (v in 1:length(varList)){
  var1    <- getVarData(inFile, varList[v])
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- varOutList[v]
  unitList                    <- c(unitList, unitVegFlux[v])
  countCol                    <- countCol + 1
}
remove(varList, varOutList, unitVegFlux)

##
inFile      <- paste0("vegfluxH2O", inFileExt)
varList     <- c("Evaporation", "Interception", "Transpiration", "GC", "GA")
unitVegFlux <- c("kg m-2 s-1", "kg m-2 s-1", "kg m-2 s-1", "m s-1", "m s-1")
for (v in 1:length(varList)){
  var1    <- getVarData(inFile, varList[v])
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- varList[v]
  unitList                    <- c(unitList, unitVegFlux[v])
  countCol                    <- countCol + 1
}
remove(varList, unitVegFlux)

var1    <- getVarData(inFile, "Runoff")
var2    <- getVarData(inFile, "Drainage")
outArr[,countCol]           <- var1 + var2
colnames(outArr)[countCol]  <- "Leaching_and_Runoff"      # kg m-2 s-1
unitList                    <- c(unitList, "kg m-2 s-1")
countCol                    <- countCol + 1

##
inFile  <- paste0("vegpoolC", inFileExt)
var1    <- getVarData(inFile, "Leaf_C")
outArr[,countCol]           <- var1
outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-2 -> kg C m-2
colnames(outArr)[countCol]  <- "Leaf_C"
unitList                    <- c(unitList, "kg C m-2")
countCol                    <- countCol + 1
leafCMol                    <- var1         # needed for calculation of foliar_N and foliar_P

var1    <- getVarData(inFile, "FineRoot_C")
outArr[,countCol]           <- var1
outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-2 -> kg C m-2
colnames(outArr)[countCol]  <- "FineRoot_C"
unitList                    <- c(unitList, "kg C m-2")
countCol                    <- countCol + 1

var1    <- getVarData(inFile, "Total_C")
outArr[,countCol]           <- var1
outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-2 -> kg C m-2
colnames(outArr)[countCol]  <- "Total_C"
unitList                    <- c(unitList, "kg C m-2")
countCol                    <- countCol + 1

##
inFile  <- paste0("vegpooldC13DC14", inFileExt)
var1    <- getVarData(inFile, "Leaf_dC13")
outArr[,countCol]           <- var1
outArr[,countCol]           <- outArr[,countCol] * molar_mass_C / 1000   # mol C m-2 -> kg C m-2
colnames(outArr)[countCol]  <- "Leaf_dC13"
unitList                    <- c(unitList, "kg C m-2")
countCol                    <- countCol + 1

##
inFile  <- paste0("vegpoolN", inFileExt)
var1    <- getVarData(inFile, "Leaf_N")
outArr[,countCol]           <- (var1 * molar_mass_N) / (leafCMol * molar_mass_C) * 48
colnames(outArr)[countCol]  <- "foliar_N"         # [%]
unitList                    <- c(unitList, "%")
countCol                    <- countCol + 1
leafNMol                    <- var1         # needed for calculation of foliar_NP and chlorophyll

##
inFile  <- paste0("vegpoolP", inFileExt)
var1    <- getVarData(inFile, "Leaf_P")
outArr[,countCol]           <- (var1 * molar_mass_P) / (leafCMol * molar_mass_C) * 48
colnames(outArr)[countCol]  <- "foliar_P"         # [%]
unitList                    <- c(unitList, "%")
countCol                    <- countCol + 1
leafPMol                    <- var1         # needed for calculation of foliar_NP
remove(leafCMol)

##
var1    <- (leafNMol * molar_mass_N) / (leafPMol * molar_mass_P)
outArr[,countCol]           <- var1
colnames(outArr)[countCol]  <- "foliar_N/P"       # [-]
unitList                    <- c(unitList, "--")
countCol                    <- countCol + 1
remove(leafPMol)

##
inFile  <- paste0("veg_canopy", inFileExt)        # canopy layers are organized in rows
cLayers <- getVarData(inFile, "CanopyLayer")      # get vector with canopy layer information for sorting data per layer

var0    <- getVarData(inFile, "f_chl")            # get data for all canopy layers
var0    <- var0 * leafNMol * 1000 * 1000 / 39.8   # chlorophyll -> by layer (1 top) chl = f_chl * leaf_nitrogen (mmol N m-2 (leaf)) *1000._wp/39.8_wp ---> µmol Chloroplast per m-2 leaf
remove(leafNMol)
for (l in 1:nCanopyLay){
  print(paste0("get and save chlorophyll for canopy layer: ", l))
  rowsLayer <- which(cLayers == l)
  var1                        <- var0[rowsLayer]
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- paste0("chlorophyll_lay", l)     # [µmol Chloroplast per m-2 leaf]
  unitList                    <- c(unitList, "µmol Chloroplast per m-2 leaf")
  countCol                    <- countCol + 1
}

var0    <- getVarData(inFile, "jmax")   # get data for all canopy layers
for (l in 1:nCanopyLay){
  print(paste0("get and save jmax for canopy layer: ", l))
  rowsLayer <- which(cLayers == l)
  var1                        <- var0[rowsLayer]
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- paste0("jmax_lay", l)     # [µmol m-2 (leaf) s-1]
  unitList                    <- c(unitList, "µmol m-2 (leaf) s-1")
  countCol                    <- countCol + 1
}

var0    <- getVarData(inFile, "vcmax")  # get data for all canopy layers
for (l in 1:nCanopyLay){
  print(paste0("get and save vcmax for canopy layer: ", l))
  rowsLayer <- which(cLayers == l)
  var1                        <- var0[rowsLayer]
  outArr[,countCol]           <- var1
  colnames(outArr)[countCol]  <- paste0("vcmax_lay", l)     # [µmol m-2 (leaf) s-1]
  unitList                    <- c(unitList, "µmol m-2 (leaf) s-1")
  countCol                    <- countCol + 1
}

###
### use only the columns of outArr with data
###

print(paste0("created an array with ", countCol-1, " columns of data."))
print(paste0("removing ", nOutCol - (countCol-1), " empty columns"))
outArr <- outArr[,1:(countCol-1)]


###
### save outArr to output file
###

print(paste0("writing output file: ", outFile))
write.table(outArr, file=outFile, append=FALSE, sep = ",", eol="\n", row.names=FALSE, col.names=TRUE)


###
### create a figure for each output variable (and layer)
###

print(paste0("creating figures across output variables, saved as: ", outFile, ".pdf"))
pdf(file=paste0(outFile, ".pdf"), paper="a4", height=8, width=8, pointsize=10, bg="white", family="Helvetica", version="1.6")
  for (c in 1:ncol(outArr)){
    plot(1:nOutRow, outArr[,c], type="p", pch=16, cex=0.6, 
         xlab="11 Years, 2009 - 2019, 365 days per year [30 min timesteps]", 
         ylab=colnames(outArr)[c], 
         main=paste0(colnames(outArr)[c], " [", unitList[c], "]"))
  }
dev.off()




## fin ##
print("Done!")
