#!/bin/bash

## copy the soil_forcing to the BGI/SLM model-input folders where the SLM_soil model expect the forcing files
## this soil_forcing needs to be created for both FLUXNET and CRUNCEPv5 using the prep_run_slm.sh with option flag_generate_soil_forcing == T
## 
## first version: Jan Engel, MPI-BGC Jena, 20170207
##


# set manually: forcing  mode with options: static tansient (see prep_run_slm)
forcingMode="static"


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
#pathServerA="$mainPathServer""FACE/"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"

# different file names for CRUNCEPv7 depending on the forcing mode
# e.g.: VAL_1901-2015.dat / VAL_1984-2013.dat (transient / static)
crunYears="_1984-2013"
if [ $forcingMode == "transient" ]
 then
  crunYears="_1901-2015"
fi

# get current date (yyyy-mm-dd) for renaming folders in user scratch
currDate=$(date +%F)


#####################################################################
## loop over potential numbers of soil layers present: 5 9 15      ##
#####################################################################
for sl in 5 9 15
 do
  echo " "
  echo " Running soil_forcing update for ""$sl"" soil layers."
  echo " "

## set the soil_forcing folder for the particular sl number at the server where the soil_forcing files are stored
#  create this folder in case it is not present
soilForce="new_soil_forcing_sl""$sl""/"
mkdir -p "$pathServerF""$soilForce"
mkdir -p "$pathServerC""$soilForce"

## set the experimentalID folder names in the user scratch/slm/ folder containing the simulation results
expNameF="exp_create_soil_forc_fluxnet_sl""$sl"
expNameCC="exp_create_soil_forc_cruncep_cue_sl""$sl"
expNameCS="exp_create_soil_forc_cruncep_spp_sl""$sl"


#################################################
## delete old and copy new soil_forcing files  ##
#################################################

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

      filenameBase=$(ls -la "$pathServerF""met_forcing/" | grep "$siteBase" | awk '{print $9}' | rev | cut -c5- |rev)
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
  echo " Removing old cruncep CUE soil_forcing files from: ""$pathServerC""$soilForce"
  echo " and copy cruncep SPP soil_forcing from :""$pathSourceBase""$expNameCS"
  for site in $(ls -d "$pathSourceBase""$expNameCS"/*/)
   do

    siteBase=$(basename $site)  # just: BBR_01 & LUE_01

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

sync
done  ## loop over soil layers

echo " "
echo "> New soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data) <"
echo " "

## fin ##
