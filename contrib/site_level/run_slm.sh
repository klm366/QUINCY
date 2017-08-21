#!/bin/bash

## 
##
## script to run simulations of SLM; i.e. one specific job (of an array, or a specific site)
##


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## arguments passed to run_slm.sh ##
jobName=$1                # a user specific alphanumeric string
simDefaultFiles=$2        # name of folder with binary namelist lctlib R-Script in $bgiScratchUserFold/slm/$jobName
singleSiteName=$3         # only provided if a single site should be simulated (not for param sensitivity), else a job array will be run


## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"
simFoldDefaultFiles="$simFold""$simDefaultFiles"   # namelist binary lctlib R-Script config-file (provided by prep_run_slm.sh)


## job array ID & job index number ##
#jobArrayID=$LSB_JOBID             # ID of the job array (assigned by OpenLava; e.g. 78435) | LSB / OpenLava internal variable
jobIndex=$LSB_JOBINDEX            # index number of the job (in the job array: 2-x)        | LSB / OpenLava internal variable


## read simulation config-file created by prep_run_slm.sh ##
siteList=$(sed -n 1p "$simFoldDefaultFiles""config.txt")    # -n  == quiet/silent
siteDataFoldForcing=$(sed -n 2p "$simFoldDefaultFiles""config.txt")   # folder at the Server containing the forcing files
slmModel=$(sed -n 3p "$simFoldDefaultFiles""config.txt")
forcingDataSet=$(sed -n 4p "$simFoldDefaultFiles""config.txt")
forcingModeFileID=$(sed -n 5p "$simFoldDefaultFiles""config.txt")
specialFluxnetOutput=$(sed -n 6p "$simFoldDefaultFiles""config.txt")  # T/F
siteDataObsFold=$(sed -n 7p "$simFoldDefaultFiles""config.txt")
useSiteSpecificAWC=$(sed -n 8p "$simFoldDefaultFiles""config.txt")    # T/F
doParamSensitivitySim=$(sed -n 9p "$simFoldDefaultFiles""config.txt" | awk '{print $1}')  # T/F
 # info: $compareSimCRUNCEPvsGFDB T/F is saved at line 8 but necessary only for the cleanup_slm_hpc_run.sh

## model forcing - define input file-names as expected by the model (see mo_atmland_forcing.f90)
forceSoilLitterFile="soil_forcing_litterfall.dat"
forceSoilPhysicFile="soil_forcing_soilphysics.dat"
forceLandMetFile="climate.dat"


# read extra config-file for parameter sensitivity simulations created by prep_run_slm.sh if applicable ##
if [ $doParamSensitivitySim == "T" ]
  then
    simParamSensiFiles=$(sed -n 9p "$simFoldDefaultFiles""config.txt" | awk '{print $2}') # name of folder containing files
    lineNumberOfSite=$(sed -n 1p "$simFold""$simParamSensiFiles""param_config.txt")
    slmParameterListDefault=$(sed -n 2p "$simFold""$simParamSensiFiles""param_config.txt")
    slmParameterSensiTable=$(sed -n 3p "$simFold""$simParamSensiFiles""param_config.txt")
fi


## further variables ##
slmNamelist="namelist.slm"
slmLctlib="lctlib_slm_nlct8.def"
slmBinary="land.x"
slmRscript="postprocessing_""$slmModel"".R"           ## just the R-Script for postprocessing output data
slmRscript4fluxfile="convert_slm2FLUXNET.R"           ## just the R-Script for converting output for FLUXNET 
slmRscript4fluxdiag="generate_FLUXNET_diagnostics.R"  ## just the R-Script for FLUXNET diagnostics 


##############################################################################
## get job / site specific data ##############################################
##############################################################################
siteData=0             # "vector" containing all data of the site as presented in the "$siteList" / all_site_list.dat
 # 1:Site-ID 2:lon 3:lat 4:PFT 5:start 6:end 7:lon_gmt 8:clay 9:silt 10:sand
 # 11:awc 12:bd 13:pH 14:taxusda 15:taxnwrb 16:lith_glim 17:LAI 18:Nleaf 19:SLA 20:Height 21:PlantYear
 # access via: "$(awk '{print $1}' <<< "$siteData")"
if [ -z "$singleSiteName" ]  # if unset (i.e. if argument $3 is not provided because a job array should be run)
  then  # normal or parameter-sensitivity simulations
    if [ $doParamSensitivitySim == T ]
      then # parameter sensitivity - a site from customSiteList
	      siteData=$(sed -n "$lineNumberOfSite"p "$siteList") # get all data of the specific site (using line number)
      else # "normal" simulations - a site from the all_site_list.dat
	      siteData=$(sed -n "$jobIndex"p "$siteList")         # get all data of the specific site (using line number)
    fi
  else  # a single site
	siteData=$(grep "$singleSiteName" "$siteList")            # get all data of the specific site using site name, e.g., [FLUX] "DE-Hai" "ZA-Kru" [CRUN] "BHS_01" "ILO_05" [FACE] "DK-AMB" "OR-ELE"
fi


##############################################################################
## job specific actions start here ###########################################
##############################################################################
# define a site-specific job-folder within the simulation folder
jobFold="$simFold""$(awk '{print $1}' <<< "$siteData")""/"
# for parameter sensitivity simulations prepend a number: 1001 == 1st set of parameter value combinations in param_table.txt
if [ $doParamSensitivitySim == "T" ]; then
  paramSimNumber=$((1000 + $jobIndex - 1))
  jobFold="$simFold""$paramSimNumber""_""$(awk '{print $1}' <<< "$siteData")""/"
fi
mkdir -p "$jobFold"
cd "$jobFold"           ## cd into job specific folder ##
jobErrFold="./error/"
mkdir -p "$jobErrFold"  # folder for error files and out files of this job

## copy files specific to this job into the job folder ##
# namelist, binary, lctlib, R-Script
cp "$simFoldDefaultFiles""$slmNamelist" ./
cp "$simFoldDefaultFiles""$slmLctlib"   ./
cp "$simFoldDefaultFiles""$slmBinary"   ./
cp "$simFoldDefaultFiles""$slmRscript"  ./
# special fluxnet-output files
if [ $specialFluxnetOutput == "T" ];then
   cp "$simFoldDefaultFiles""$slmRscript4fluxfile"  ./
   cp "$simFoldDefaultFiles""$slmRscript4fluxdiag"  ./
   cp "$siteDataObsFold""$(awk '{print $1}' <<< "$siteData")"*.nc ./obs.nc
fi
# forcing file(s)
forcingSiteFileID="$(awk '{print $1}' <<< $siteData)"
if [ $forcingDataSet == "cruncep" ]
 then
  forcingSiteFileID=$(cut -c1-3 <<< $forcingSiteFileID) 
fi
case "$slmModel" in
"land")
  cp "$siteDataFoldForcing""$forcingSiteFileID""$forcingModeFileID"*.dat      ./"$forceLandMetFile"
  ;;
"plant")
  cp "$siteDataFoldForcing""$forcingSiteFileID""$forcingModeFileID"*.dat      ./"$forceLandMetFile"
  ;;
"canopy")
  cp "$siteDataFoldForcing""$forcingSiteFileID""$forcingModeFileID"*.dat      ./"$forceLandMetFile"
  ;;
"soil")
  cp "$siteDataFoldForcing""$forcingSiteFileID""$forcingModeFileID"*_lit.dat      ./"$forceSoilLitterFile"
  cp "$siteDataFoldForcing""$forcingSiteFileID""$forcingModeFileID"*_sph.dat      ./"$forceSoilPhysicFile"
  ;;
esac
# parameter sensitivity files
if [ $doParamSensitivitySim == "T" ]; then
  slmParameterList="param.list"
  cp ../"$simParamSensiFiles""$slmParameterListDefault"           ./"$slmParameterList"
fi

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
    echo "EXIT!  no such PFT is known by the run_slm.sh script: ""$pftName"
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

# 1:Site-ID 2:lon 3:lat 4:PFT 5:start 6:end 7:lon_gmt 8:clay 9:silt 10:sand 11:awc 12:bd 13:pH 14:taxusda 15:taxnwrb 16:lith_glim 17:LAI 18:Nleaf 19:SLA 20:Height 21:PlantYear
changeNamelistValue longitude             "$(awk '{print $2}'  <<< "$siteData")"
changeNamelistValue latitude              "$(awk '{print $3}'  <<< "$siteData")"
changeNamelistValue plant_functional_type_id "$pftID"
changeNamelistValue soil_clay             "$(awk '{print $8}'  <<< "$siteData")"
changeNamelistValue soil_silt             "$(awk '{print $9}'  <<< "$siteData")"
changeNamelistValue soil_sand             "$(awk '{print $10}' <<< "$siteData")"
changeNamelistValue bulk_density          "$(awk '{print $12}' <<< "$siteData")"
changeNamelistValue soil_ph               "$(awk '{print $13}' <<< "$siteData")"
changeNamelistValue usda_taxonomy_class   "$(awk '{print $14}' <<< "$siteData")"
changeNamelistValue nwrb_taxonomy_class   "$(awk '{print $15}' <<< "$siteData")"
changeNamelistValue glim_lithology_class  "$(awk '{print $16}' <<< "$siteData")"
[ $useSiteSpecificAWC == "T" ] && \
   changeNamelistValue soil_awc_prescribe "$(awk '{print $11}' <<< "$siteData")"   # this is equal to: if [ ... ]; then ...; fi
# special fluxnet output #
if [ $specialFluxnetOutput == "T" ]
  then
    year1=$(awk '{print $5}' <<< "$siteData")
    year2=$(awk '{print $6}' <<< "$siteData")
    forcingPeriodLength=$(($year2-$year1+1))
    changeNamelistValue fluxnet_output_length "$forcingPeriodLength"
fi

## manipulate the param.list if parameter sensitivity simulations should be done ##
if [ $doParamSensitivitySim == "T" ]; then

  # define function that sets the actual values per parameter in ./param.list
  function changeParamListValue {
    sed -i /"$1"/s/ACTIVE/"$2"/ ./"$slmParameterList"   # replace the placeholder "ACTIVE" for parameter $1 with value $2
  return 0                                            # get return value with: return_value=$?
  }
  # number of parameters
  paraNumber=$(sed -n 1p ../"$simParamSensiFiles""$slmParameterSensiTable" | grep -o " " | wc -l) # count the number of separators between columns

  for ((i = 0 ; i < $paraNumber ; i++))
    do
      paraCount=$(($i +2))
      paraName=$(sed -n 1p ../"$simParamSensiFiles""$slmParameterSensiTable" | awk -v c="$paraCount" '{print $c}' )
      paraValue=$(sed -n "$jobIndex"p ../"$simParamSensiFiles""$slmParameterSensiTable" | awk -v c="$paraCount" '{print $c}' ) 
      changeParamListValue $paraName $paraValue
  done
fi


## sync all data with the filesystem ##
sync

## run the simulation ##
./"$slmBinary"
sync

## run the postprocessing R Script ##
Rscript ./"$slmRscript"
sync

if [ $specialFluxnetOutput == "T" ]; then 
   Rscript ./"$slmRscript4fluxfile" $(awk '{print $2}' <<< "$siteData") $(awk '{print $3}' <<< "$siteData") $(awk '{print $5}' <<< "$siteData")
   Rscript ./"$slmRscript4fluxdiag" $(awk '{print $5}' <<< "$siteData") $(awk '{print $6}' <<< "$siteData")
fi

## cleanup some job specific files ##
rm ./*.dat  # forcing file(s)
rm ./"$slmRscript"
rm ./"$slmBinary"
[ $specialFluxnetOutput == "T" ] && rm ./"$slmRscript4fluxfile" ./"$slmRscript4fluxdiag"  # equal to: if [ ... ]; then ...; fi
sync

## sort job specific files ##
mkdir -p ./input
mkdir -p ./output
mkdir -p ./postprocessing
sync
mv ./"$slmNamelist"         ./input
mv ./"$slmLctlib"           ./input
if [ $doParamSensitivitySim == "T" ]; then
  mv ./"$slmParameterList"  ./input
fi
if [ $specialFluxnetOutput == "T" ]; then 
  mv ./summary*             ./postprocessing
  mv ./statistics*          ./postprocessing
  mv ./*.nc                 ./output 
fi
mv ./*.txt                  ./output 
if [ "$slmModel" == "land" ] || [ "$slmModel" == "plant" ] || [ "$slmModel" == "canopy" ] # || [ "$slmModel" == "soil" ]
 then
  mv ./*.pdf                  ./postprocessing
fi

## sync all written data with file system ##
sync

## end run_slm.sh ##
