#!/bin/bash

## 20160601
##
## script to run simulations of SLM; i.e. one specific job (of an array, or a specific site)
##


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## arguments passed to run_slm.sh ##
jobName=$1                        # a user specific alphanumeric string
simDefaultFiles=$2                # name of folder with binary namelist R-Script in $bgiScratchUserFold/$jobName
singleSiteName=$3                 # only provided if a single site should be simulated, else a job array will be run


## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"
simFoldDefaultFiles="$simFold""$simDefaultFiles"   # namelist binary R-Script


## job array ID & job index number ##
#jobArrayID=$LSB_JOBID             # ID of the job array (assigned by OpenLava; e.g. 78435) | LSB / OpenLava internal variable
jobIndex=$LSB_JOBINDEX            # index number of the job (in the job array: 2-x)        | LSB / OpenLava internal variable


## read simulation config-file created by prep_run_slm.sh ##
siteList=$(sed -n 1p "$simFoldDefaultFiles""config.txt")    # -n  == quiet/silent
siteDataFold=$(sed -n 2p "$simFoldDefaultFiles""config.txt")
slmModel=$(sed -n 3p "$simFoldDefaultFiles""config.txt")
specialFluxnetOutput=$(sed -n 4p "$simFoldDefaultFiles""config.txt")
siteDataObsFold=$(sed -n 5p "$simFoldDefaultFiles""config.txt")
useSiteSpecificAWC=$(sed -n 6p "$simFoldDefaultFiles""config.txt")

## further variables ##
slmNamelist="namelist.slm"
slmBinary="land.x"
slmRscript="postprocessing_""$slmModel"".R"           ## just the R-Script for postprocessing output data
slmRscript4fluxfile="convert_slm2FLUXNET.R"           ## just the R-Script for converting output for FLUXNET 
slmRscript4fluxdiag="generate_FLUXNET_diagnostics.R"  ## just the R-Script for FLUXNET diagnostics 


##############################################################################
## get job / site specific data ##############################################
##############################################################################
siteData=0             # containing all data of the site as presented in the "$siteList" / all_site_list.dat
# 1:Site-ID 2:lon 3:lat 4:PFT 5:start 6:end 7:lon_gmt 8:clay 9:silt 10:sand 11:awc 12:pH 13:LAI 14:Nleaf 15:SLA 16:Height 17:PlantYear
# access via: "$(awk '{print $1}' <<< "$siteData")"
if [ -z "$singleSiteName" ]  # if unset (i.e. if $3 is not provided because a job array should be run)
  then  # one site of a job array
	siteData=$(sed -n "$jobIndex"p "$siteList")       # get all data of the specific site (using line number)
  else  # a single site
	siteData=$(grep "$singleSiteName" "$siteList")          # get all data of the specific site (using site name, e.g. DE-Bay)
fi


##############################################################################
## job specific actions start here ###########################################
##############################################################################
jobFold="$simFold""$(awk '{print $1}' <<< "$siteData")""/"
mkdir -p "$jobFold"
cd "$jobFold"
jobErrFold="./error/"
mkdir -p "$jobErrFold"  # folder for error files and out files of this job

## copy files specific to this job into the job folder ##
# namelist binary R-Script all_site_list.dat
cp "$simFoldDefaultFiles""$slmNamelist" ./
cp "$simFoldDefaultFiles""$slmBinary"   ./
cp "$simFoldDefaultFiles""$slmRscript"  ./
if [ $specialFluxnetOutput == "yes" ];then
   cp "$simFoldDefaultFiles""$slmRscript4fluxfile" ./
   cp "$simFoldDefaultFiles""$slmRscript4fluxdiag"  ./
   cp "$siteDataObsFold""$(awk '{print $1}' <<< "$siteData")"*.nc ./obs.nc
fi
cp "$siteDataFold""$(awk '{print $1}' <<< "$siteData")"*.dat ./climate.dat

## convert PFT name into ID ##
pftName="$(awk '{print $4}' <<< "$siteData")"
pftID=0
case "$pftName" in
  TrBE)
    pftID=1 ;;
  TeBE)
    pftID=2 ;;
  TrBR)
    pftID=3 ;;
  TeBS)
    pftID=4 ;;
  BBS)
    pftID=4 ;;
  BNE)
    pftID=5 ;;
  TeNE)
    pftID=5 ;;
  BNS)
    pftID=6 ;;
  TeH)
    pftID=7 ;;
  TeC)
    pftID=7 ;;
  TrH)
    pftID=8 ;;
  TrC)
    pftID=8 ;;
  *)
    echo "EXIT!  no such PFT is known by the run_slm.sh script"
    exit 
esac


## manipulate the namelist - parameters specific to this job / site ##
 # function: pass the variable_name ($1) & the new_value ($2)
 # $slmNamelist refers to the namelist file name, assuming the namelist file is located in the current folder !
function changeNamelistValue {
oldValue=$(grep -w $1 ./"$slmNamelist" | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 ./"$slmNamelist" | cut -d : -f 1)          # get line number of variable $1
sed -i "$lineNum"s/$oldValue/$2/ ./"$slmNamelist"                # replace the old value in line $lineNum with $2
return 0                                                         # get return value with: return_value=$?
}
 # examples:
 #changeNamelistValue simulation_length_number 300
 #changeNamelistValue output_interval_flux monthly    # quotation marks not needed
 #changeNamelistValue flag_calc_soil_moisture .FALSE.
 #changeNamelistValue soil_theta_prescribe 122.7

changeNamelistValue longitude "$(awk '{print $2}' <<< "$siteData")"
changeNamelistValue latitude  "$(awk '{print $3}' <<< "$siteData")"
changeNamelistValue plant_functional_type_id "$pftID"
[ $useSiteSpecificAWC == "yes" ] && \
   changeNamelistValue soil_awc_prescribe "$(awk '{print $11}' <<< "$siteData")"   # this is equal to: if [ ... ]; then ...; fi
# special fluxnet output #
if [ $specialFluxnetOutput == "yes" ]
  then
    year1=$(awk '{print $5}' <<< "$siteData")
    year2=$(awk '{print $6}' <<< "$siteData")
    forcingPeriodLength=$(($year2-$year1+1))
    changeNamelistValue fluxnet_output_length "$forcingPeriodLength"
fi
sync


## run the simulation ##
./"$slmBinary"

## sync all written data with file system ##
sync

## run the R Script ##
Rscript ./"$slmRscript"
sync

if [ $specialFluxnetOutput == "yes" ]; then 
   Rscript ./"$slmRscript4fluxfile" $(awk '{print $2}' <<< "$siteData") $(awk '{print $3}' <<< "$siteData") $(awk '{print $5}' <<< "$siteData")
   Rscript ./"$slmRscript4fluxdiag" $(awk '{print $5}' <<< "$siteData") $(awk '{print $6}' <<< "$siteData")
fi

## cleanup some job specific files ##
rm ./climate.dat
rm ./"$slmRscript"
rm ./"$slmBinary"
[ $specialFluxnetOutput == "yes" ] && rm ./"$slmRscript4fluxfile" ./"$slmRscript4fluxdiag"  # equal to: if [ ... ]; then ...; fi

## sort job specific files ##
mkdir -p ./input
mkdir -p ./output
mkdir -p ./postprocessing
mv ./"$slmNamelist"  ./input
if [ $specialFluxnetOutput == "yes" ]; then 
  mv ./summary*      ./postprocessing
  mv ./statistics*   ./postprocessing
  mv ./*.nc ./output 
fi
mv ./*.txt           ./output 
mv ./*.pdf           ./postprocessing

## sync all written data with file system ##
sync

## end run_slm.sh ##
