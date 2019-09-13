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

## copy the just-created soil_forcing (usecase 01) to the particular folders in the sites sets (BGI/data/Model/SLM1/input/point/)
##  considering soil forcing for 15 soil layers only
##
## the soil forcing mey be created for FLUXNET and CRUNCEP cue & spp and FACE and SPP_1685:
##  using the prep_run_slm.sh with the usecase 01, which sets default values and the option flag_generate_soil_forcing == T
## 



# no need to set any options, just run this script to copy new soil forcing available



########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## get username as needed for the location of the simulation results
userName="$USER"

## paths to the folders at the server
mainPathServer="/Net/Groups/BGI/data/Model/SLM1/input/point/"
pathServerF="$mainPathServer""FLUXNET2/"
pathServerC="$mainPathServer""CRUNCEPv7/"
pathServerA="$mainPathServer""FACE/"
pathServerS="$mainPathServer""SPP_1685/"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"

# met-forcing years included into the file names for CRUNCEP
crunYears="_1984-2013"

# get current date (yyyy-mm-dd) for renaming folders in user scratch after soil forcing was copied
currDate=$(date +%F)


#################################################
## define main variables and create folders    ##
#################################################

## set number of soil layers (15 is the default; currently we use only 15 sl soil-forcing)
sl=15
echo " "
echo " Running soil_forcing update for ""$sl"" soil layers for all site-sets for which new soil forcing is available."
echo " "

## set the soil_forcing folder for the particular sl number at the server where the soil_forcing files are stored
#  create this folder in case it is not present, does preserve folders and symlinks to folders present
soilForce="soil_forcing_sl""$sl""/"
mkdir -p "$pathServerF""$soilForce"
mkdir -p "$pathServerC""$soilForce"
mkdir -p "$pathServerA""$soilForce"
mkdir -p "$pathServerS""$soilForce"

## set the experimentalID folder names in the user scratch/slm/ folder containing the simulation results
expNameF="exp_create_soil_forc_fluxnet_sl""$sl"
expNameCC="exp_create_soil_forc_cruncep_cue_sl""$sl"
expNameCS="exp_create_soil_forc_cruncep_spp_sl""$sl"
expNameA="exp_create_soil_forc_face_sl""$sl"
expNameS="exp_create_soil_forc_spp1685_sl""$sl"


########################################################################################
## delete old and copy new soil_forcing files, only if new soil forcing was created   ##
########################################################################################

#############################
## Fluxnet                 ##
#############################
if [ -d "$pathSourceBase""$expNameF" ] # if folder is present, i.e., new soil_forcing files were created
 then

  echo " Removing old Fluxnet soil_forcing from: ""$pathServerF""$soilForce"
  rm -f "$pathServerF""$soilForce"*lit.dat
  rm -f "$pathServerF""$soilForce"*sph.dat
  sync && sleep 0.5

  echo " Copying Fluxnet soil_forcing from :""$pathSourceBase""$expNameF"
  for site in $(ls -d "$pathSourceBase""$expNameF"/*/)
   do

    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...

    if [ $siteBase != "x_binary" ]
    then

      # grep the met_forcing filename with e.g. DE-Hai, select only the first line (static forcing, DE-Hai_s_2000-2006) and remove "_s"
      filenameBase=$(ls -la "$pathServerF""met_forcing/" | grep "$siteBase" | head -n1 | awk '{print $9}' | rev | cut -c5- | rev | sed "s/_s//")
      echo "  site: "$filenameBase

      # copy new soil_forcing files
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerF""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerF""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameF"  "$pathSourceBase""$expNameF"_"$currDate"
fi     # if new soil_forcing was created


#############################
## Cruncep CUE             ##
#############################
if [ -d "$pathSourceBase""$expNameCC" ] # if folder is present, i.e., new soil_forcing files were created
 then

  echo " "
  echo " Removing old cruncep CUE soil_forcing files from: ""$pathServerC""$soilForce"
  echo " and copy cruncep CUE soil_forcing from :""$pathSourceBase""$expNameCC"
  for site in $(ls -d "$pathSourceBase""$expNameCC"/*/)
   do

    siteBase=$(basename $site)  # just: FLA_01, NAU_01, NAU_02, TUM_01, ...

    if [ $siteBase != "x_binary" ]
    then
      # only for the _01 sites
      numStrg=$(echo $siteBase | rev | cut -c1-2 | rev)  
      if [ "$numStrg" == "01" ]
       then

        # get a short version of siteBase FLA_01 -> FLA (which is consistent with the name of the met forcing file)
        siteBaseShort=$(echo $siteBase | cut -c1-3)
        echo  "  site: ""$siteBaseShort""_""$numStrg"

        # delete old soil_forcing files of this site
        rm -f "$pathServerC""$soilForce""$siteBaseShort"*lit.dat
        rm -f "$pathServerC""$soilForce""$siteBaseShort"*sph.dat
        sync && sleep 0.2

        # copy new soil_forcing files
        cp "$pathSourceBase""$expNameCC""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_lit.dat"
        cp "$pathSourceBase""$expNameCC""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_sph.dat"

      fi # only _01 sites
    fi   # not x_binary
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameCC" "$pathSourceBase""$expNameCC"_"$currDate"
fi     # if new soil_forcing was created


#############################
## Cruncep SPP             ##
#############################
if [ -d "$pathSourceBase""$expNameCS" ] # if folder is present, i.e., new soil_forcing files were created
 then

  echo " "
  echo " Removing old cruncep SPP soil_forcing files from: ""$pathServerC""$soilForce"
  echo " and copy cruncep SPP soil_forcing from :""$pathSourceBase""$expNameCS"
  for site in $(ls -d "$pathSourceBase""$expNameCS"/*/)
   do

    siteBase=$(basename $site)  # just: BBR_01 & LUE_01 & VES_01

    if [ $siteBase != "x_binary" ]
    then
        # the two SPP sites have a *_01 name, so do not need to sort out the *_01 folder as needed for CUE
        # get a short version of siteBase BBR_01 -> BBR (which is the name of the forcing file)
        siteBaseShort=$(echo $siteBase | cut -c1-3)
        echo  "  site: ""$siteBaseShort""_01"

        # delete old soil_forcing files of this site
        rm -f "$pathServerC""$soilForce""$siteBaseShort"*lit.dat
        rm -f "$pathServerC""$soilForce""$siteBaseShort"*sph.dat
        sync && sleep 0.2

        # echo "  site: "$filenameBase
        cp "$pathSourceBase""$expNameCS""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_lit.dat"
        cp "$pathSourceBase""$expNameCS""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_sph.dat"

    fi # not x_binary
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameCS" "$pathSourceBase""$expNameCS"_"$currDate"
fi     # if new soil_forcing was created


#############################
## Face                    ##
#############################
if [ -d "$pathSourceBase""$expNameA" ] # if folder is present, i.e., new soil_forcing files were created
 then

  echo " Removing old Face soil_forcing from: ""$pathServerA""$soilForce"
  rm -f "$pathServerA""$soilForce"*lit.dat
  rm -f "$pathServerA""$soilForce"*sph.dat
  sync && sleep 0.5

  echo " Copying Face soil_forcing from :""$pathSourceBase""$expNameA"
  # loop over all sites for which soil_forcing was created
  for site in $(ls -d "$pathSourceBase""$expNameA"/*/)
   do

    siteBase=$(basename $site)  # just: DK-AMB, DK-ELE, OR-AMB and OR-ELE

    if [ $siteBase != "x_binary" ]
    then

      # get the name of the soil-forcing files from the met-forcing files
      # in met_forcing folder two files for each site are present, one for each: static and transient forcing
      #   the 2nd is used, which is the static one, with a later date (e.g. DK-AMB_1901-2015 & DK-AMB_1996-2007)
      filenameBase=$(ls -la "$pathServerA""met_forcing/" | grep "$siteBase" | tail -n1 | awk '{print $9}' | rev | cut -c5- | rev)
      echo "  site: "$filenameBase

      # copy new soil_forcing files
      cp "$pathSourceBase""$expNameA""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerA""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameA""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerA""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameA"  "$pathSourceBase""$expNameA"_"$currDate"
fi     # if new soil_forcing was created


#############################
## SPP_1685                ##
#############################
if [ -d "$pathSourceBase""$expNameS" ] # if folder is present, i.e., new soil_forcing files were created
 then

  echo " Removing old SPP_1685 soil_forcing from: ""$pathServerS""$soilForce"
  rm -f "$pathServerS""$soilForce"*lit.dat
  rm -f "$pathServerS""$soilForce"*sph.dat
  sync && sleep 0.5

  echo " Copying SPP_1685 soil_forcing from :""$pathSourceBase""$expNameS"
  # loop over all sites for which soil_forcing was created
  for site in $(ls -d "$pathSourceBase""$expNameS"/*/)
   do

    siteBase=$(basename $site)  # just: VES, LUE, BBR, ... 

    if [ $siteBase != "x_binary" ]
    then

      # get the name of the soil-forcing files from the met-forcing files
      # in met_forcing folder two files for each site are present, one for each: static and transient forcing
      #   the 2nd is used, which is the static one, with a later date (e.g. BBR_1901-2015.dat & BBR_1984-2013.dat)
      filenameBase=$(ls -la "$pathServerS""met_forcing/" | grep "$siteBase" | tail -n1 | awk '{print $9}' | rev | cut -c5- | rev)
      echo "  site: "$filenameBase

      # copy new soil_forcing files
      cp "$pathSourceBase""$expNameS""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerS""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameS""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerS""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameS"  "$pathSourceBase""$expNameS"_"$currDate"
fi     # if new soil_forcing was created


sync

echo " "
echo "> New soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data) <"
echo " "

## fin ##
