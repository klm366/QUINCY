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
## works only with FLUXNET2/3 site-sets
## 

# no need to set any options, just run this script to copy new transient soil forcing available

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
pathServerF3="$mainPathServer""FLUXNET3/"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"

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
mkdir -p "$pathServerF3""$soilForce"

## set the experimentalID folder names in the user scratch/slm/ folder containing the simulation results
expNameF="exp_create_soil_forc_trans_fluxnet_sl""$sl"
expNameF3="exp_create_soil_forc_trans_fluxnet3_sl""$sl"

########################################################################################
## delete old and copy new soil_forcing files, only if new soil forcing was created   ##
########################################################################################

#############################
## Fluxnet 2               ##
#############################
if [ -d "$pathSourceBase""$expNameF" ]; then # if folder is present, i.e., new transient soil_forcing files were created
  # removing old transient forcing files
  echo " Removing old Fluxnet transient soil_forcing from: ""$pathServerF""$soilForce"
  rm -f "$pathServerF""$soilForce"*_t_*lit.dat
  rm -f "$pathServerF""$soilForce"*_t_*sph.dat
  sync && sleep 0.5
  # copy & rename new transient forcing files
  echo " Copying Fluxnet transient soil_forcing from :""$pathSourceBase""$expNameF"
  for site in $(ls -d "$pathSourceBase""$expNameF"/*/)
   do
    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...
    if [ $siteBase != "x_binary" ]; then
      # grep the static met_forcing filename from e.g. "AT-Neu_s_", replace _s_ by _t_ and remove .dat => AT-Neu_t_2002-2012
      # this provides the correct forcing period written out by u08, i.e., the last full forcing period in 1901-2015
      filenameBase=$(ls -la "$pathServerF""met_forcing/" | grep "$siteBase""_s_" | awk '{print $9}' | sed "s/_s_/_t_/" | sed "s/.dat//")
      echo "  site: "$filenameBase
      # copy new soil_forcing files
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerF""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerF""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameF"  "$pathSourceBase""$expNameF"_"$currDate"
fi  # if new soil_forcing had been created

#############################
## Fluxnet 3               ##
#############################
if [ -d "$pathSourceBase""$expNameF3" ]; then # if folder is present, i.e., new transient soil_forcing files were created
  # removing old transient forcing files
  echo " Removing old Fluxnet3 transient soil_forcing from: ""$pathServerF3""$soilForce"
  rm -f "$pathServerF3""$soilForce"*_t_*lit.dat
  rm -f "$pathServerF3""$soilForce"*_t_*sph.dat
  sync && sleep 0.5
  # copy & rename new transient forcing files
  echo " Copying Fluxnet3 transient soil_forcing from :""$pathSourceBase""$expNameF3"
  for site in $(ls -d "$pathSourceBase""$expNameF3"/*/)
   do
    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...
    if [ $siteBase != "x_binary" ]; then
      # grep the static met_forcing filename from e.g. "AT-Neu_s_", replace _s_ by _t_ and remove .dat => AT-Neu_t_2002-2012
      # this provides the correct forcing period written out by u08, i.e., the last full forcing period in 1901-2015
      filenameBase=$(ls -la "$pathServerF3""met_forcing/" | grep "$siteBase""_s_" | awk '{print $9}' | sed "s/_s_/_t_/" | sed "s/.dat//")
      echo "  site: "$filenameBase
      # copy new soil_forcing files
      cp "$pathSourceBase""$expNameF3""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerF3""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameF3""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerF3""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done # loop over all sites for which soil_forcing was created
  ## appending current date to simulation results folders at your scratch
  mv "$pathSourceBase""$expNameF3"  "$pathSourceBase""$expNameF3"_"$currDate"
fi  # if new soil_forcing had been created

######################################
## done message                     ##
######################################
echo " "
echo "> New transient fluxnet/3 soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data) <"
echo " "

## fin ##
