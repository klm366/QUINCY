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

## 
## script to run simulations of SLM; i.e. this script is called individually for each specific job (of an array, or a specific site)
## 
## this script is copied to /Net/Groups/BGI/scratch/""$username""/slm/$jobName/default_files/ by the prep_run_slm.sh
##  and the bsup runs it from there
##


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## arguments passed to run_slm.sh ##
jobName=$1                # a user specific alphanumeric string: "exp_""$experimentalID"
simDefaultFiles=$2        # name of folder with binary namelist lctlib R-Script in $bgiScratchUserFold/slm/$jobName
singleSiteName=$3         # only provided if a single site should be simulated (not for param sensitivity), else a job array will be run

## set a random sleep time (1 to 90 sec) to not let each instance of this run_slm.sh run at exactly the same time
# this helps to reduce the max load of canossa.bgc-jena.mpg.de thereby improving BGC-Filesystem stability
wait=$(( $RANDOM % 90 )) 
sleep $wait

## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"
simFoldDefaultFiles="$simFold""$simDefaultFiles"  # this runscript, namelist, binary, lctlib, R-Script for postprocessing, config-file (provided by prep_run_slm.sh)


## job array ID & job index number ##
#jobArrayID=$LSB_JOBID             # ID of the job array (assigned by OpenLava; e.g. 78435) | LSB / OpenLava internal variable
jobIndex=$LSB_JOBINDEX            # index number of the job (in the job array: 2-x)        | LSB / OpenLava internal variable


## read simulation config-file created by prep_run_slm.sh ##
# sed -n  == quiet/silent
siteList=$(sed -n 1p "$simFoldDefaultFiles""config.txt")  
siteDataFoldForcing=$(sed -n 2p "$simFoldDefaultFiles""config.txt")   # folder at the Server containing the forcing files
slmModel=$(sed -n 3p "$simFoldDefaultFiles""config.txt")
forcingDataSet=$(sed -n 4p "$simFoldDefaultFiles""config.txt")
forcingModeFileID=$(sed -n 5p "$simFoldDefaultFiles""config.txt")
specialFluxnetOutput=$(sed -n 6p "$simFoldDefaultFiles""config.txt")                                  # T/F
siteDataObsFold=$(sed -n 7p "$simFoldDefaultFiles""config.txt")
useSiteSpecificAWC=$(sed -n 8p "$simFoldDefaultFiles""config.txt")                                    # T/F
doParamSensitivitySim=$(sed -n 9p "$simFoldDefaultFiles""config.txt" | awk '{print $1}')              # T/F
# info: $compareSimCRUNCEPvsGFDB T/F is saved at line 10 but not needed here; only for the cleanup_slm_hpc_run.sh
fluxnetTransientSim=$(sed -n 11p "$simFoldDefaultFiles""config.txt" | awk '{print $1}')               # T/F
transient_simulation_start_year=$(sed -n 11p "$simFoldDefaultFiles""config.txt" | awk '{print $2}')   # year at which transient simulation starts (default 1901)
createSoilForcFlag=$(sed -n 12p "$simFoldDefaultFiles""config.txt")                                   # T/F

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

## define a function that checks whether folder is present replacing 'sync'
# function argument $1 is the folder to be checked
function checkDirPresent {
 if [ ! -d "$1" ]; then
   sleep 2
   if [ ! -d "$1" ]; then
     sleep 2
     mkdir -p "$1"
     if [ ! -d "$1" ]; then
       sleep 2
       sync
       if [ ! -d "$1" ]; then
         echo " EXIT! Directory ""$1"" could not be created (run_slm.sh)."
       fi
     fi
   fi
 fi
return 0                                                         # get return value with: return_value=$?
}

## define a site-specific job-folder within the simulation folder
siteName="$(awk '{print $1}' <<< "$siteData")"
jobFold="$simFold""$siteName""/"
# for parameter sensitivity simulations prepend a number (where 1001 == 1st set of parameter value combinations in parameter_table.txt)
if [ $doParamSensitivitySim == "T" ]; then
  paramSimNumber=$((1000 + $jobIndex - 1))
  jobFold="$simFold""$paramSimNumber""_""$(awk '{print $1}' <<< "$siteData")""/"
fi

## create job specific folder using the absolute path ##
# "$bgiScratchUserFold""$jobName""/""$(awk '{print $1}' <<< "$siteData")""/"
mkdir -p "$jobFold"
checkDirPresent "$jobFold"

## cd into job specific folder ##
cd "$jobFold"           ## cd into job specific folder ##
jobErrFold="./error/"
jobOutFold="./std_out/"
mkdir -p "$jobErrFold"  # folder for error files of this job
mkdir -p "$jobOutFold"  # folder for files with std output of this job
checkDirPresent "$jobErrFold"
checkDirPresent "$jobOutFold"

## copy files specific to this job into the job folder ##
# namelist, binary, lctlib, R-Script
cp "$simFoldDefaultFiles""$slmNamelist" ./
cp "$simFoldDefaultFiles""$slmLctlib"   ./
cp "$simFoldDefaultFiles""$slmBinary"   ./
cp "$simFoldDefaultFiles""$slmRscript"  ./
# special fluxnet-output files (only if not used for soil forcing creation)
if [ $specialFluxnetOutput == "T" ] && [ "$createSoilForcFlag" == "F" ] && [ "$slmModel" == "land" ]; then
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
# copy the "parameter sensitivity namelist" i.e. the input file for parameter values for SLM as used in mo_shared_set_parameters
if [ $doParamSensitivitySim == "T" ]; then
  slmParameterList="parameter_slm_run.list"       # this name must be identical to the filename set in mo_shared_set_parameters:t_parameters_from_file
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
sleep 2s && sync


######################################################################
## manipulate the namelist - parameters specific to this job / site ##
######################################################################
 # function: pass the variable_name ($1) & the new_value ($2)
 # $slmNamelist refers to the namelist file name, assuming the namelist file is located in the current folder !
function changeNamelistValue {
oldValue=$(grep -w $1 ./"$slmNamelist" | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 ./"$slmNamelist" | cut -d : -f 1)          # get line number of variable $1
sed -i "$lineNum"s/$oldValue/$2/ ./"$slmNamelist"                # replace the old value in line $lineNum with $2
# explicitly testing and error catching to avoid using sync
sleep 0.1
newValue=$(grep -w $1 ./"$slmNamelist" | awk -F "=" '{print $2}' | tr -d '"') # get the new value of this variable
if [ "$2" != "$newValue" ]; then
  sleep 2
  newValue=$(grep -w $1 ./"$slmNamelist" | awk -F "=" '{print $2}' | tr -d '"') # get the new value of this variable
  if [ "$2" != "$newValue" ]; then
    sync
  fi
fi
return 0                                                         # get return value with: return_value=$?
}
 # examples:
 #changeNamelistValue simulation_length_number 300
 #changeNamelistValue output_interval_flux monthly    # quotation marks not needed
 #changeNamelistValue flag_calc_soil_moisture .FALSE.
 #changeNamelistValue soil_theta_prescribe 122.7

# 1:Site-ID 2:lon 3:lat 4:PFT 5:start 6:end 7:lon_gmt 8:clay 9:silt 10:sand 
# 11:awc 12:bd 13:pH 14:taxusda 15:taxnwrb 16:lith_glim 17:LAI 18:Nleaf 19:SLA 20:Height 
# 21:PlantYear 22:soil_p_depth 23:soil_p_labile 24:soil_p_slow 25:soil_p_occluded 26:soil_p_primary 27:qmax_org_fine_particle 
changeNamelistValue longitude               "$(awk '{print $2}'  <<< "$siteData")"
changeNamelistValue latitude                "$(awk '{print $3}'  <<< "$siteData")"
changeNamelistValue plant_functional_type_id "$pftID"
changeNamelistValue soil_clay               "$(awk '{print $8}'  <<< "$siteData")"
changeNamelistValue soil_silt               "$(awk '{print $9}'  <<< "$siteData")"
changeNamelistValue soil_sand               "$(awk '{print $10}' <<< "$siteData")"
changeNamelistValue bulk_density            "$(awk '{print $12}' <<< "$siteData")"
changeNamelistValue soil_ph                 "$(awk '{print $13}' <<< "$siteData")"
changeNamelistValue usda_taxonomy_class     "$(awk '{print $14}' <<< "$siteData")"
changeNamelistValue nwrb_taxonomy_class     "$(awk '{print $15}' <<< "$siteData")"
changeNamelistValue lithology_class         "$(awk '{print $16}' <<< "$siteData")"
changeNamelistValue lai_site_specific       "$(awk '{print $17}' <<< "$siteData")"
changeNamelistValue stand_replacing_year    "$(awk '{print $21}' <<< "$siteData")"
changeNamelistValue soil_p_depth            "$(awk '{print $22}' <<< "$siteData")"
changeNamelistValue soil_p_labile           "$(awk '{print $23}' <<< "$siteData")"
changeNamelistValue soil_p_slow             "$(awk '{print $24}' <<< "$siteData")"
changeNamelistValue soil_p_occluded         "$(awk '{print $25}' <<< "$siteData")"
changeNamelistValue soil_p_primary          "$(awk '{print $26}' <<< "$siteData")"
changeNamelistValue qmax_org_fine_particle  "$(awk '{print $27}' <<< "$siteData")"
# use site-specific AWC #
if [ $useSiteSpecificAWC == "T" ]; then
   changeNamelistValue soil_awc_prescribe "$(awk '{print $11}' <<< "$siteData")" 
fi
# special fluxnet output #
# calc forcing period length for namelist - needed also for soil forcing creation
if [ $specialFluxnetOutput == "T" ]; then
  # startyear
  year1=$(awk '{print $5}' <<< "$siteData")
  # in case of transient simulations: get startyear from "site-info file for static simulations"
  if [ $fluxnetTransientSim == "T" ]; then
    if [ $forcingDataSet == "fluxnet" ]; then
      year1=$(grep "$forcingSiteFileID" "/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat" | awk '{print $5}')
    elif [ $forcingDataSet == "fluxnet3" ]; then
      year1=$(grep "$forcingSiteFileID" "/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET3/flux3_all_sites_1991-2014_list.dat" | awk '{print $5}')
    fi
  fi
  # lastyear
  year2=$(awk '{print $6}' <<< "$siteData")
  # period
  forcingPeriodLength=$(($year2-$year1+1))
  changeNamelistValue fluxnet_output_length "$forcingPeriodLength"
fi
# transient simulations with fluxnet sites #
if [ $fluxnetTransientSim == "T" ]; then
  lastDataYearSiteSpecific="$(awk '{print $6}' <<< "$siteData")"
  simulationLengthSiteSpecific=$((lastDataYearSiteSpecific - transient_simulation_start_year + 1))
  changeNamelistValue simulation_length_number "$simulationLengthSiteSpecific"
  changeNamelistValue output_end_last_day_year "$simulationLengthSiteSpecific"
fi


## manipulate the parameter_slm.list if parameter sensitivity simulations should be done ##
if [ $doParamSensitivitySim == "T" ]; then

  # define function that sets the actual values per parameter in ./parameter_slm.list
  function changeParamListValue {
    lineNum=$(grep -wn $1 "$slmParameterList" | cut -d : -f 1)  # get line number of $paramName
    sed -i "$lineNum"s/ACTIVE/"$2"/ ./"$slmParameterList"       # replace the placeholder "ACTIVE" for the parameter in line $lineNum with value $2
    sleep 0.4
  return 0                                                      # get return value with: return_value=$?
  }
  # number of parameters
  paraNumber=$(sed -n 1p ../"$simParamSensiFiles""$slmParameterSensiTable" | grep -o " " | wc -l) # count the number of separators between columns

  # loop over parameter names in the header of $slmParameterSensiTable
  for ((i = 0 ; i < $paraNumber ; i++))
    do
      # loop from 0 to nParameter-1; the 1st parameter is in column 2
      paraCount=$(($i +2))
      # get the parameter name
      paraName=$(sed -n 1p ../"$simParamSensiFiles""$slmParameterSensiTable" | awk -v c="$paraCount" '{print $c}' )
      # get the value of this parameter for this specific run
      paraValue=$(sed -n "$jobIndex"p ../"$simParamSensiFiles""$slmParameterSensiTable" | awk -v c="$paraCount" '{print $c}' ) 
      # set the value for this parameter (replacing the key word 'ACTIVE')
      changeParamListValue $paraName $paraValue
  done
fi

## sync all data with the filesystem ##
sleep 5s && sync


#############################################################################################
## move the job dir to the local /scratch at the HPC node and work from there              ##
#############################################################################################
## test how much space is left on the local /scratch of the node
scratchSpace=$(df /scratch/ | grep "scratch" | awk '{print $4}')
# get the name of the node
nodeName=$(uname -n)
echo "simulation runs at node: ""$nodeName"
# assume a QUINCY job might take max 4 GB space (for an u82 job 3 GB are used), calc the overall min space needed
# default is 32 GB (8 jobs at all queue nodes) [1 GB = 1 B * 1028 * 1028]
minSpace="34000000"
# 512 GB at access io1 io2 io3 where 128 jobs can run
if [ "$nodeName" == "access" ] || [ "$nodeName" == "io1" ] || [ "$nodeName" == "io2" ] || [ "$nodeName" == "io3" ]; then
  minSpace="544000000"
# 192 GB at io4-6 where 48 jobs can run
elif [ "$nodeName" == "io4" ] || [ "$nodeName" == "io5" ] || [ "$nodeName" == "io6" ]; then
  minSpace="204000000"
# 128 GB at master where 32 jobs can run
elif [ "$nodeName" == "master" ]; then
  minSpace="136000000"
# 80 GB at iof-ioj where 20 jobs can run
elif [ "$nodeName" == "iof" ] || [ "$nodeName" == "iog" ] || [ "$nodeName" == "ioh" ] || [ "$nodeName" == "ioi" ] || [ "$nodeName" == "ioj" ]; then
  minSpace="85000000"
# 64 GB at io7-ioe where 16 jobs can run
elif [ "$nodeName" == "io7" ] || [ "$nodeName" == "io8" ] || [ "$nodeName" == "io9" ] || [ "$nodeName" == "ioa" ] || [ "$nodeName" == "iob" ] || [ "$nodeName" == "ioc" ] || [ "$nodeName" == "iod" ] || [ "$nodeName" == "ioe" ]; then
  minSpace="68000000"
fi

localScratch="FALSE"
if [ "$scratchSpace" -gt "$minSpace" ]
 then
  localScratch="TRUE"
  echo "running the simulation at the local /scratch of the node"
  echo " "

  ## working in the scratch folder at the node
  mkdir -p /scratch/"$username"
  mkdir -p /scratch/"$username"/"$jobName""$jobIndex""$siteName"
  rm -rf   /scratch/"$username"/"$jobName""$jobIndex""$siteName"/*     # remove any files just in the unlikely case the folder did already exist 
  cd       /scratch/"$username"/"$jobName""$jobIndex""$siteName"

  ## copy all files (but not the folder) here
  cp -r   "$jobFold"* ./ 
  ## make sure all folders and files are copied to local scratch
  sleep 10s && sync
  ## remove all copied files from network/scratch
  rm -rf  "$jobFold"*
fi


## get start-time of simulation
startTime=$(date -u +%F_%X_%Z)

## run the simulation ##
./"$slmBinary"

## get finish-time of simulation
finishTime=$(date -u +%F_%X_%Z)
sleep 2

## run the postprocessing R Script ##
Rscript ./"$slmRscript"
sync

if [ $specialFluxnetOutput == "T" ] && [ "$createSoilForcFlag" == "F" ] && [ "$slmModel" == "land" ]; then 
   Rscript ./"$slmRscript4fluxfile" $(awk '{print $2}' <<< "$siteData") $(awk '{print $3}' <<< "$siteData") $(awk '{print $5}' <<< "$siteData")
   Rscript ./"$slmRscript4fluxdiag" $(awk '{print $5}' <<< "$siteData") $(awk '{print $6}' <<< "$siteData")
fi

## cleanup some job specific files ##
rm ./*.dat                  #  forcing file(s) (climate.dat for SLM land)
rm ./"$slmRscript"          #  R script for postprocessing
rm ./"$slmBinary"           #  the binary
if [ $specialFluxnetOutput == "T" ] && [ "$createSoilForcFlag" == "F" ] && [ "$slmModel" == "land" ]; then 
  rm ./"$slmRscript4fluxfile"
  rm ./"$slmRscript4fluxdiag"
fi
sleep 2

## sort job specific files ##
mkdir -p ./input
mkdir -p ./output
mkdir -p ./postprocessing
checkDirPresent ./input
checkDirPresent ./output
checkDirPresent ./postprocessing

mv ./"$slmNamelist"         ./input               #  namelist.slm
mv ./"$slmLctlib"           ./input               #  usually this would be: lctlib_slm_nlct8.def
mv ./*.txt                  ./output              #  the output files SLM created at runtime
mv ./*.pdf                  ./postprocessing      #  for all SLM models postprocessing is available & executed
if [ $doParamSensitivitySim == "T" ]; then
  mv ./"$slmParameterList"  ./input
fi
if [ $specialFluxnetOutput == "T" ] && [ "$createSoilForcFlag" == "F" ] && [ "$slmModel" == "land" ]; then 
  mv ./output/summary*      ./postprocessing
  mv ./output/statistics*   ./postprocessing
  mv ./*.nc                 ./output 
fi

## write start & finish times of this job to file
echo "#####"                              >  "$jobOutFold"start_finish_sim_time.txt
echo "$jobIndex"                          >> "$jobOutFold"start_finish_sim_time.txt
echo $(awk '{print $1}' <<< "$siteData")  >> "$jobOutFold"start_finish_sim_time.txt
echo "start:  ""$startTime"               >> "$jobOutFold"start_finish_sim_time.txt
echo "finish: ""$finishTime"              >> "$jobOutFold"start_finish_sim_time.txt
echo "  "                                 >> "$jobOutFold"start_finish_sim_time.txt

## sync all written data with file system ##
sync

## cd back into the network/scratch & copy data there
## if 
if [ "$localScratch" == "TRUE" ]
 then
  cd "$jobFold"
  cp -r   /scratch/"$username"/"$jobName""$jobIndex""$siteName"/* ./
  sync         # sync the copied data, so that it is actually written to the harddrive
  sleep 2      # an extra layer of security to ensure the data has been copied before removing job-folder from local scratch
  # remove job-folder from local scratch
  rm -rf  /scratch/"$username"/"$jobName""$jobIndex""$siteName"
 else
  echo " Did not run the simulation at the local /scratch of the node because there was not enough space available!"
  echo " "
fi

## end run_slm.sh ##
