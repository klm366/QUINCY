#!/bin/bash

## copy the soil_forcing to the BGI/SLM model-input folders where the SLM_soil model expect the forcing files
## this soil_forcing needs to be created for both FLUXNET and CRUNCEPv5 using the prep_run_slm.sh with option flag_generate_soil_forcing == T
## 
## first version: Jan Engel, MPI-BGC Jena, 20170207
##


# CRUNCEPv7 forcing  mode (see prep_run_slm)
# options: static tansient
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
soilForce="soil_forcing/"

## source folders - where the simulation results have been stored
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"
pathSourceBase="/Net/Groups/BGI/scratch/""$userName""/slm/"
expNameF="exp_create_soil_forc_fluxnet"
expNameCC="exp_create_soil_forc_cruncep_cue"
expNameCS="exp_create_soil_forc_cruncep_spp"

# different file names for CRUNCEPv7 depending on the forcing mode
# e.g.: VAL_1901-2015.dat / VAL_1984-2013.dat (transient / static)
crunYears="_1984-2013"
if [ $forcingMode == "transient" ]
 then
  crunYears="_1901-2015"
fi


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

    siteBase=$(basename $site)  # just: DE-Hai, IT-Amp, US-Syv, ...

    if [ $siteBase != "x_binary" ]
    then

      #echo $siteBase
      filenameBase=$(ls -la "$pathServerF""met_forcing/" | grep "$siteBase" | awk '{print $9}' | rev | cut -c5- |rev)
      #echo $filenameBase

      echo "  site: "$filenameBase
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerF""$soilForce""$filenameBase""_lit.dat"
      cp "$pathSourceBase""$expNameF""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerF""$soilForce""$filenameBase""_sph.dat"
    fi
    sync
  done

## Cruncep CUE
echo " Copying Cruncep soil_forcing from :""$pathSourceBase""$expNameCC"
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
        # #echo $siteBase
        # filenameBase=$(ls -la "$pathServerC""met_forcing/" | grep "$siteBaseShort" | awk '{print $9}' | rev | cut -c5- |rev)
        # #echo $filenameBase

        # echo "  site: "$filenameBase
        echo  "  site: ""$siteBaseShort""_""$numStrg"
        cp "$pathSourceBase""$expNameCC""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_lit.dat"
        cp "$pathSourceBase""$expNameCC""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_sph.dat"

      fi # only _01 sites
    fi # not x_binary
    sync
  done

## Cruncep SPP
echo " Copying Cruncep soil_forcing from :""$pathSourceBase""$expNameCS"
  for site in $(ls -d "$pathSourceBase""$expNameCS"/*/)
   do

    siteBase=$(basename $site)  # just: BBR_01 & LUE_01

    if [ $siteBase != "x_binary" ]
    then
      # # only for the _01 sites
      # numStrg=$(echo $siteBase | rev | cut -c1-2 | rev)  
      # if [ "$numStrg" == "01" ]
      #  then

        # get a short version of siteBase FLA_01 -> FLA (which is the name of the forcing file)
        siteBaseShort=$(echo $siteBase | cut -c1-3)
        # #echo $siteBase
        # filenameBase=$(ls -la "$pathServerC""met_forcing/" | grep "$siteBaseShort" | awk '{print $9}' | rev | cut -c5- |rev)
        # #echo $filenameBase

        # echo "  site: "$filenameBase
        echo  "  site: ""$siteBaseShort""_01"
        cp "$pathSourceBase""$expNameCS""/""$siteBase"/output/veglitterfall_timestep.txt "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_lit.dat"
        cp "$pathSourceBase""$expNameCS""/""$siteBase"/output/soilphysics_timestep.txt   "$pathServerC""$soilForce""$siteBaseShort""$crunYears""_sph.dat"

      # fi # only _01 sites
    fi # not x_binary
    sync
  done


## appending current date to simulation results folders at your scratch
currDate=$(date +%F)
mv "$pathSourceBase""$expNameF"  "$pathSourceBase""$expNameF""$currDate"
mv "$pathSourceBase""$expNameCC" "$pathSourceBase""$expNameCC""$currDate"
mv "$pathSourceBase""$expNameCS" "$pathSourceBase""$expNameCS""$currDate"

sync
echo "> New soil_forcing was copied to the ""$mainPathServer""* soil_forcing folders (replacing old data) <"


## fin ##
