#!/bin/bash

## copy the soil_forcing to the BGI/SLM model-input folders where the SLM_soil model expect the forcing files
## this soil_forcing needs to be created for both FLUXNET and CRUNCEPv5 using the prep_run_slm.sh with option flag_generate_soil_forcing == T
## 
## first version: Jan Engel, MPI-BGC Jena, 20170207
##


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## get username as needed for the location of the simulation results
userName="$USER"

## paths to the folders at the server
pathServerF="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET/"
pathServerC="/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv5/"
soilForce="soil_forcing/"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"
expNameF="exp_create_soil_forc_f"
expNameC="exp_create_soil_forc_c"

## remove old soil_forcing
rm -f "$pathServerF""$soilForce"*lit.dat
rm -f "$pathServerF""$soilForce"*sph.dat
rm -f "$pathServerC""$soilForce"*lit.dat
rm -f "$pathServerC""$soilForce"*sph.dat
sync
echo " Old soil_forcing files deleted"


## copy new soil_forcing files ##

## Fluxnet
echo " Copying Fluxnet soil_forcing from :""$pathSourceBase""$expNameF"
  for site in $(ls -d "$pathSourceBase""$expNameF"/*/)
   do

    siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67

    if [ $siteBase != "x_binary" ]
    then

      #echo $siteBase
      filenameBase=$(ls -la "$pathServerF" | grep "$siteBase" | awk '{print $9}' | rev | cut -c5- |rev)
      #echo $filenameBase

      echo "  site: "$filenameBase
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/veglitterfall_tile1_timestep.txt "$pathServerF""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/soilphysics_tile1_timestep.txt   "$pathServerF""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done

## Cruncep
echo " Copying Cruncep soil_forcing from :""$pathSourceBase""$expNameC"
  for site in $(ls -d "$pathSourceBase""$expNameC"/*/)
   do

    siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67

    if [ $siteBase != "x_binary" ]
    then

      #echo $siteBase
      filenameBase=$(ls -la "$pathServerC" | grep "$siteBase" | awk '{print $9}' | rev | cut -c5- |rev)
      #echo $filenameBase

      echo "  site: "$filenameBase
      cp "$pathSourceBase""$expNameC""/""$siteBase"/output/veglitterfall_tile1_timestep.txt "$pathServerC""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameC""/""$siteBase"/output/soilphysics_tile1_timestep.txt   "$pathServerC""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done

## appending current date to simulation results folders at your scratch
currDate=$(date +%F)
mv "$pathSourceBase""$expNameF" "$pathSourceBase""$expNameF""$currDate"
mv "$pathSourceBase""$expNameC" "$pathSourceBase""$expNameC""$currDate"

sync
echo "> New soil_forcing was copied to the BGI/SLM model-input folders (replacing old data) <"


## fin ##
