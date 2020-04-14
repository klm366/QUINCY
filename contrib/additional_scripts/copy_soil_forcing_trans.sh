#!/bin/bash

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

## copy the just-created transient soil_forcing (usecase 08) to the particular folders in the sites sets (BGI/data/Model/SLM1/input/point/)
##  some more detailed docu is presented in the usecase 08 file 08_generate_flux_trans_soil_forcing_def.incl
##  considering soil forcing for 15 soil layers only
##
## works only with FLUXNET2/3 site-sets (also u08 does work only with FLUXNET2/3)
## 

# no need to set any options, just run this script to copy new transient soil forcing available

########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################

## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## message
echo " "
echo " Copying new transient soil_forcing (15 soil layers) for all site-sets for which new soil forcing is available."
echo " "

# get current date (yyyy-mm-dd) for renaming folders in user scratch after soil forcing was copied
currDate=$(date +%F)

## get username as needed for the location of the simulation results
userName="$USER"

## model name dir
modelNameSourceDir="slm"
modelNameTargetDir="SLM1"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/""$modelNameSourceDir""/"

## define the directory names without revision number (rX) in scratch/user/"$modelNameSourceDir"/ containing the simulation results
#  revision number is detected per site-set
expNameF="exp_create_soil_forc_trans_fluxnet_sl15_r"
expNameF3="exp_create_soil_forc_trans_fluxnet3_sl15_r"

## paths to the directories the soil-forcing may be copied to per site-set
mainPathServer="/Net/Groups/BGI/data/Model/""$modelNameTargetDir""/input/point/"
pathServerF="$mainPathServer""FLUXNET2/"
pathServerF3="$mainPathServer""FLUXNET3/"

## set the base-name of the soil_forcing directory valid for all site-sets (differing for CRUNCEP cue/spp)
soilForceDirBase="soil_forcing_trans_sl15_r"

#################################################################################################
## copy new soil_forcing files (delete old of same rev.) only if new soil forcing was created  ##
#################################################################################################

#############################
## FLUXNET2                ##
#############################
siteSet="FLUXNET2"
expName="$expNameF"
pathServer="$pathServerF"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are two or more directories of ""$expName"" - remove the old one and re-run this script."
    exit 1
  fi
  # get revision number the soil_forcing was created with (only the number, i.e., all behind the "sl15_r")
  # complete the dir name of the soil-forcing created in scratch/user/"$modelNameSourceDir"/
  revNum=$(ls -lA "$pathSourceBase" | grep "$expName" | awk -F'_sl15_r' '{print $2}')
  expName="$expName""$revNum"
  # compile path and dir name for the soil-forcing to be copied to
  sfTargetPathDir="$pathServer""$soilForceDirBase""$revNum""/"
  # dir with met-forcing used to create this soil-forcing
  metForcPathDir="$pathServer""met_forcing/"
  
  # delete previously created soil-forcing with this rev number, if present
  if [ -d "$sfTargetPathDir" ]; then
    echo " Removing old ""$siteSet"" soil_forcing: ""$sfTargetPathDir"
    rm -rf "$sfTargetPathDir"
    sleep 0.5 && sync
  fi
  # create new dir for new soil-forcing
  mkdir -p "$sfTargetPathDir"
  # create an info dir within and copy of the 1st site namelist & lctlib, and the postprocessing_all_sites.pdf file
  mkdir -p "$sfTargetPathDir""info"
  firstSiteDir=$(ls -d "$pathSourceBase""$expName"/*/ | head -n 1)
  firstSiteDir=$(basename "$firstSiteDir")
  cp "$pathSourceBase""$expName""/""$firstSiteDir""/input/"*  "$sfTargetPathDir""info""/"
  if [ -e "$pathSourceBase""$expName""/postprocessing_all_sites.pdf" ]; then
    cp "$pathSourceBase""$expName""/postprocessing_all_sites.pdf"  "$sfTargetPathDir""info""/"
  fi

  # copy & rename soil-forcing files
  echo " Copying ""$siteSet"" soil_forcing from : ""$pathSourceBase""$expName"
  echo "  to: ""$sfTargetPathDir"
  for site in $(ls -d "$pathSourceBase""$expName"/*/)
   do
    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...
    if [ $siteBase != "x_binary" ]; then
      # get the basic soil-forcing filename from static! met-forcing filename with start- & end-year; DE-Hai_s_2000-2006.dat -> DE-Hai_2000-2006
      # note: transient soil-forcing was created using transient met-forcing, but covers only the static met-forcing period
      filenameBase=$(ls -lA "$metForcPathDir" | grep "$siteBase""_s_" | awk '{print $9}' | sed "s/_s//" | sed "s/.dat//")
      echo "  site: "$filenameBase
      # copy new soil_forcing files
      cp "$pathSourceBase""$expName""/""$siteBase"/output/veglitterfall_timestep.txt "$sfTargetPathDir""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expName""/""$siteBase"/output/soilphysics_timestep.txt   "$sfTargetPathDir""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at scratch/user/"$modelNameSourceDir"/
  mv "$pathSourceBase""$expName"  "$pathSourceBase""copied_""$expName"_"$currDate"
else
  echo " NO new soil-forcing for: ""$siteSet"
fi     # if new soil_forcing was created

#############################
## FLUXNET3                ##
#############################
siteSet="FLUXNET3"
expName="$expNameF3"
pathServer="$pathServerF3"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are two or more directories of ""$expName"" - remove the old one and re-run this script."
    exit 1
  fi
  # get revision number the soil_forcing was created with (only the number, i.e., all behind the "sl15_r")
  # complete the dir name of the soil-forcing created in scratch/user/"$modelNameSourceDir"/
  revNum=$(ls -lA "$pathSourceBase" | grep "$expName" | awk -F'_sl15_r' '{print $2}')
  expName="$expName""$revNum"
  # compile path and dir name for the soil-forcing to be copied to
  sfTargetPathDir="$pathServer""$soilForceDirBase""$revNum""/"
  # dir with met-forcing used to create this soil-forcing
  metForcPathDir="$pathServer""met_forcing/"
  
  # delete previously created soil-forcing with this rev number, if present
  if [ -d "$sfTargetPathDir" ]; then
    echo " Removing old ""$siteSet"" soil_forcing: ""$sfTargetPathDir"
    rm -rf "$sfTargetPathDir"
    sleep 0.5 && sync
  fi
  # create new dir for new soil-forcing
  mkdir -p "$sfTargetPathDir"
  # create an info dir within and copy of the 1st site namelist & lctlib, and the postprocessing_all_sites.pdf file
  mkdir -p "$sfTargetPathDir""info"
  firstSiteDir=$(ls -d "$pathSourceBase""$expName"/*/ | head -n 1)
  firstSiteDir=$(basename "$firstSiteDir")
  cp "$pathSourceBase""$expName""/""$firstSiteDir""/input/"*  "$sfTargetPathDir""info""/"
  if [ -e "$pathSourceBase""$expName""/postprocessing_all_sites.pdf" ]; then
    cp "$pathSourceBase""$expName""/postprocessing_all_sites.pdf"  "$sfTargetPathDir""info""/"
  fi

  # copy & rename soil-forcing files
  echo " Copying ""$siteSet"" soil_forcing from : ""$pathSourceBase""$expName"
  echo "  to: ""$sfTargetPathDir"
  for site in $(ls -d "$pathSourceBase""$expName"/*/)
   do
    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...
    if [ $siteBase != "x_binary" ]; then
      # get the basic soil-forcing filename from static! met-forcing filename with start- & end-year; DE-Hai_s_2000-2006.dat -> DE-Hai_2000-2006
      # note: transient soil-forcing was created using transient met-forcing, but covers only the static met-forcing period
      filenameBase=$(ls -lA "$metForcPathDir" | grep "$siteBase""_s_" | awk '{print $9}' | sed "s/_s//" | sed "s/.dat//")
      echo "  site: "$filenameBase
      # copy new soil_forcing files
      cp "$pathSourceBase""$expName""/""$siteBase"/output/veglitterfall_timestep.txt "$sfTargetPathDir""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expName""/""$siteBase"/output/soilphysics_timestep.txt   "$sfTargetPathDir""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at scratch/user/"$modelNameSourceDir"/
  mv "$pathSourceBase""$expName"  "$pathSourceBase""copied_""$expName"_"$currDate"
else
  echo " NO new soil-forcing for: ""$siteSet"
fi     # if new soil_forcing was created


########################################################################################
## final sync and message                                                             ##
########################################################################################
sync
echo " "
echo "> New transient soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data with same revision number) <"
echo " "

## fin ##
