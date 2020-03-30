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

## copy the just-created static soil_forcing (usecase 01) to the particular folders in the sites sets (BGI/data/Model/SLM1/input/point/)
##  some more detailed docu is presented in the usecase 01 file 01_generate_stat_soil_forcing_def.incl
##  covering soil forcing for 15 soil layers only
##
## the soil forcing may be created for FLUXNET2/3 and CRUNCEP cue & spp and FACE and SPP_1685:
##  using the prep_run_slm.sh with the usecase 01, which set default values and the option flag_generate_soil_forcing == T
## 

# no need to set any options, just run this script to copy the new static soil-forcing available

########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################

## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## message
echo " "
echo " Running soil_forcing update for 15 soil layers for all site-sets for which new soil forcing is available."
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
expNameF="exp_create_soil_forc_stat_fluxnet_sl15_r"
expNameF3="exp_create_soil_forc_stat_fluxnet3_sl15_r"
expNameA="exp_create_soil_forc_stat_face_sl15_r"
expNameS="exp_create_soil_forc_stat_spp1685_sl15_r"
expNameCC="exp_create_soil_forc_stat_cruncep_cue_sl15_r"
expNameCS="exp_create_soil_forc_stat_cruncep_spp_sl15_r"

## paths to the directories the soil-forcing may be copied to per site-set
mainPathServer="/Net/Groups/BGI/data/Model/""$modelNameTargetDir""/input/point/"
pathServerF="$mainPathServer""FLUXNET2/"
pathServerA="$mainPathServer""FACE/"
pathServerS="$mainPathServer""SPP_1685/"
pathServerF3="$mainPathServer""FLUXNET3/"
pathServerC="$mainPathServer""CRUNCEPv7/"

## set the base-name of the soil_forcing directory valid for all site-sets, differing for CRUNCEP cue/spp
soilForceDirBase="soil_forcing_stat_sl15_r"
soilForceCruncepCueDirBase="soil_forcing_stat_cue_sl15_r"
soilForceCruncepSppDirBase="soil_forcing_stat_spp_sl15_r"

# met-forcing years included into the soil-forcing file names for CRUNCEP cue/spp
sfCruncepYears="_1984-2013"

########################################################################################
## delete old and copy new soil_forcing files, only if new soil forcing was created   ##
########################################################################################

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
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
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
      # get the basic soil-forcing filename from met-forcing filename with start- & end-year; DE-Hai_s_2000-2006.dat -> DE-Hai_2000-2006
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
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
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
      # get the basic soil-forcing filename from met-forcing filename with start- & end-year; DE-Hai_s_2000-2006.dat -> DE-Hai_2000-2006
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
## FACE                    ##
#############################
siteSet="FACE"
expName="$expNameA"
pathServer="$pathServerA"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
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
    siteBase=$(basename $site)  # just: DK-AMB, DK-ELE, OR-AMB and OR-ELE
    if [ $siteBase != "x_binary" ]; then
      # get the basic soil-forcing filename from met-forcing filename with start- & end-year; DK-AMB_1996-2007.dat -> DK-AMB_1996-2007
      #   the 2nd file is used, which is the static one with a later start-year (e.g. DK-AMB_1901-2015 & DK-AMB_1996-2007)
      filenameBase=$(ls -lA "$metForcPathDir" | grep "$siteBase" | tail -n1 | awk '{print $9}' | sed "s/.dat//")
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
## SPP_1685                ##
#############################
siteSet="SPP_1685"
expName="$expNameS"
pathServer="$pathServerS"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
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
    siteBase=$(basename $site)  # just: VES, LUE, BBR, CON, MTF
    if [ $siteBase != "x_binary" ]; then
      # get the basic soil-forcing filename from met-forcing filename with start- & end-year; BBR_1984-2013.dat -> BBR_1984-2013
      #   the 2nd file is used, which is the static one with a later start-year (e.g. BBR_1901-2015.dat & BBR_1984-2013.dat)
      filenameBase=$(ls -lA "$metForcPathDir" | grep "$siteBase" | tail -n1 | awk '{print $9}' | sed "s/.dat//")
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
## CRUNCEP CUE             ##
#############################
siteSet="CRUNCEP_CUE"
expName="$expNameCC"
pathServer="$pathServerC"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
    exit 1
  fi
  # get revision number the soil_forcing was created with (only the number, i.e., all behind the "sl15_r")
  # complete the dir name of the soil-forcing created in scratch/user/"$modelNameSourceDir"/
  revNum=$(ls -lA "$pathSourceBase" | grep "$expName" | awk -F'_sl15_r' '{print $2}')
  expName="$expName""$revNum"
  # compile path and dir name for the soil-forcing to be copied to
  sfTargetPathDir="$pathServer""$soilForceCruncepCueDirBase""$revNum""/"
  # no met_forcing dir needed here because the site directories all get the same start- & end-year: "$sfCruncepYears"
  
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
  # ! only _01 sites !
  echo " Copying ""$siteSet"" soil_forcing of the _01 sites from : ""$pathSourceBase""$expName"
  echo "  to: ""$sfTargetPathDir"
  for site in $(ls -d "$pathSourceBase""$expName"/*/)
   do
    siteBase=$(basename $site)  # just: FLA_01, NAU_01, NAU_02, TUM_01, ...
    if [ "$siteBase" != "x_binary" ]; then
      siteNumStrg=$(echo "$siteBase" | rev | cut -c1-2 | rev)  
      if [ "$siteNumStrg" == "01" ]; then
        siteBaseShort=$(sed "s/_01//" <<< "$siteBase")
        echo "  site: ""$siteBase"
        # copy new soil_forcing files
        cp "$pathSourceBase""$expName""/""$siteBase"/output/veglitterfall_timestep.txt "$sfTargetPathDir""$siteBaseShort""_""$sfCruncepYears""_lit.dat"
        cp "$pathSourceBase""$expName""/""$siteBase"/output/soilphysics_timestep.txt   "$sfTargetPathDir""$siteBaseShort""_""$sfCruncepYears""_sph.dat"
      fi
      sync
    fi
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at scratch/user/"$modelNameSourceDir"/
  mv "$pathSourceBase""$expName"  "$pathSourceBase""copied_""$expName"_"$currDate"
else
  echo " NO new soil-forcing for: ""$siteSet"
fi     # if new soil_forcing was created

#############################
## CRUNCEP SPP             ##
#############################
siteSet="CRUNCEP_SPP"
expName="$expNameCS"
pathServer="$pathServerC"
# if directory (with any revision number attached) is present, i.e., new stat soil_forcing was created
if [ -d "$pathSourceBase""$expName"* ]; then
  # test  if there is only one "$expName"* directory present (does not work with > 1)
  numExpDir=$(ls -lA "$pathSourceBase" | grep "$expName" | wc -l)
  if [ "$numExpDir" -gt 1 ]; then
    echo " EXIT! There are more than one directories of ""$expName"" - remove the old one and re-run this script."
    exit 1
  fi
  # get revision number the soil_forcing was created with (only the number, i.e., all behind the "sl15_r")
  # complete the dir name of the soil-forcing created in scratch/user/"$modelNameSourceDir"/
  revNum=$(ls -lA "$pathSourceBase" | grep "$expName" | awk -F'_sl15_r' '{print $2}')
  expName="$expName""$revNum"
  # compile path and dir name for the soil-forcing to be copied to
  sfTargetPathDir="$pathServer""$soilForceCruncepSppDirBase""$revNum""/"
  # no met_forcing dir needed here because the site directories all get the same start- & end-year: "$sfCruncepYears"
  
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
  # ! only _01 sites !
  echo " Copying ""$siteSet"" soil_forcing of the _01 sites from : ""$pathSourceBase""$expName"
  echo "  to: ""$sfTargetPathDir"
  for site in $(ls -d "$pathSourceBase""$expName"/*/)
   do
    siteBase=$(basename $site)  # just: BBR_01 & LUE_01 & VES_01 & CON_01 & MTF_01
    if [ "$siteBase" != "x_binary" ]; then
      siteNumStrg=$(echo "$siteBase" | rev | cut -c1-2 | rev)  
      if [ "$siteNumStrg" == "01" ]; then
        siteBaseShort=$(sed "s/_01//" <<< "$siteBase")
        echo "  site: ""$siteBase"
        # copy new soil_forcing files
        cp "$pathSourceBase""$expName""/""$siteBase"/output/veglitterfall_timestep.txt "$sfTargetPathDir""$siteBaseShort""_""$sfCruncepYears""_lit.dat"
        cp "$pathSourceBase""$expName""/""$siteBase"/output/soilphysics_timestep.txt   "$sfTargetPathDir""$siteBaseShort""_""$sfCruncepYears""_sph.dat"
      fi
      sync
    fi
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
echo "> New soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data) <"
echo " "

## fin ##
