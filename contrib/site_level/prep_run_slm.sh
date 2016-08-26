#!/bin/bash

## 20160811
##
## > prepare and start run_slm.sh script to run SLM simulations across multiple sites
## > the prep_run_slm.sh & run_slm-sh should be located in trunk/contrib/site_level/ (same for any branch)
## > you can run this script from any folder, i.e. you have not to cd into contrib/site_level/ 
##   (this is because this script knows about its location, due to the command: 'cd "$( dirname "$0" )"')
## 

########################################################################################################
####### simulation specific options you may want to change #############################################
########################################################################################################

# folder name where data will be stored is: /Net/Groups/BGI/scratch/USER/slm/exp_'experimentalID'
experimentalID=""         # set something useful for you (allowed only: alphanumeric '_' '-'); if empty ( experimentalID="" ) current date is used
# select sites: FLUXNET or CRUNCEPv5
siteDataSet="f"          # f: FLUXNET or c: CRUNCEPv5  (use small caps!)
# want a custom site list? -> default name is: "my_site_list.dat"; must be located in same folder as this script (usually "contrib/site_level/")
customSiteList="no"          # set either to yes or no (small caps)
customSiteListFileName="my_fluxnet_site_list.dat"    # you can specify a different file name other than my_site_list.dat; leave blank to use the default file name
# specify SLM model -> the parameter 'slm_model_name' in the namelist.slm is changed accordingly (so, no need to manually change it in the namelist)
slmModel="land"             # land, plant are currently available
# run only a single site? then provide a site name, else write: singleSiteName=""
singleSiteName=""            # one site name, e.g. "DE-Bay"; must be consistent with the parameter value: siteDataSet
# run Fluxnet simulations with data output at each timestep for the last complete forcing period (see land_model_driver.f90)
specialFluxnetOutput="no"    # set either to yes or no (small caps)

## namelist.slm options ##
# flags: T (True) or F (False) [large caps]
# plant-soil model options
flag_calc_soil_uptake="T"
flag_calc_soil_moisture="T"
flag_base_soil_model="T"
# canopy model options
flag_t_jmax_acclimation="F"
flag_optimal_Nfraction="F"
# plant growth model options
flag_optimal_Nleaf="F"
flag_optimal_allocation="F"
# vegetation dynamics options
flag_tile_veg_dynamics="F"
flag_start_bareground="F"

# simulation length & output
simulation_length_number=30     # 1 to ...500 ? [simulation_length_unit defaults to "y" == years]
output_interval_pool="weekly"   # timestep daily weekly monthly yearly
output_interval_flux="weekly"   # timestep daily weekly monthly yearly
output_start_first_day_year=1   # 1 to ... [lower than or equal with output_end_last_day_year ]
output_end_last_day_year=30     # equal with or larger than output_start_first_day_year

# forcing input
soil_awc_prescribe=300.0        # float !
useSiteSpecificAWC="no"      # yes/no to use AWC as in site-list file (soil_awc_prescribe is used otherwise)
soil_theta_prescribe=100.0      # float !

## details about the HPC, if unsure stick with these default values ##
queue="BTM"                  # available queues: all io BTM mpi_large
emailNotification="no"       # set either to yes or no (small caps)
memory=1000000               # memory per job you may need; provided in Kb (1 Gb == 1000000 Kb)
numberCpuMax=80              # maximum number of simulations started in parallel

########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## prepare Validation of user-defined model-options ##
# Function to validate model-options provided by the user
# $1 must be a part of $2 list, else the script is stopped
# example: validateUserConfigValues "$queue" "$validateQueue"
function validateUserConfigValues {
validate="$(grep -w "$1" <<< "$2")"
if [ -z "${validate// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  echo "EXIT! > ""$1"" <  is not a valid value for parameter ""$3"", please change and restart prep_run_slm.sh"
  exit
 else
  if [ -z $1 ]
   then
    echo "EXIT!  The value for parameter ""$3"" is empty but must contain a valid string, please change and restart prep_run_slm.sh"
    exit
  fi
fi
return 0
}
# lists to validate against
validatesiteDataSet="f c"
validateCustomSiteList="yes no"
validateSlmModel="land plant"  # soil canopy
validateQueue="all io BTM mpi_large"
validateEmailNotification="yes no"
validateSpecialFluxnetOutput="yes no"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validateOptionForcing="constant generate readfile"
validateUseSiteSpecificAWC="yes no"

## set and create basic folders and files ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
mkdir -p "$bgiScratchUserFold"
slmBinary="../../x86_64-intel/bin/land.x"                  ## just one binary for all models
slmNamelistFilename="namelist.slm"
slmNamelist="../namelist/""$slmNamelistFilename"           ## just one namelist for all models
slmRscript="../r_script/postprocessing_""$slmModel"".R"    ## just the R-Script for postprocessing output data
slmRscript4fluxfile="../r_script/convert_slm2FLUXNET.R"    ## just the R-Script for converting output for FLUXNET 
slmRscript4fluxdiag="../r_script/generate_FLUXNET_diagnostics.R"    ## just the R-Script for FLUXNET diagnostics 

## experimental ID & simulation specific folders at scratch ##
# if experimental ID is empty, use date as: YYYYMMDD_HHMMSS
if [ -z $experimentalID ]
 then
  experimentalID=$(date +%Y%m%d_%H%M%S)
fi
# define job name
jobName="exp_""$experimentalID"  # name of the job(array) -> must start with a character make work (see below): bsub -w "ended($jobName)"
simFold="$bgiScratchUserFold""$jobName""/"
simDefaultFiles="default_files/"
simFoldDefaultFiles="$simFold""$simDefaultFiles"
simFoldOut="$simFold""out""/"
simFoldErr="$simFold""error""/"
# check if folder with this experimentalID already exists
if [ -d "$simFold" ]
 then
  echo "Note: you already used th experimentalID "$jobName"! Should previous data be deleted?"
  echo " Type either y (yes delete), a (to add to this data), or n (no, stop this script), and hit Enter."
  echo " Default is n."
  read proceed
   if [ -z "$proceed" ] # if $proceed is empty -> default
    then
     proceed="n"
   fi
  case "$proceed" in
    n)
      echo "Ok, the script has stopped!"
      exit
      ;;
    y)
      echo "Ok, all previous data within the folder with this experimentalID will be deleted, and the script processes as normal."
      rm -r "$simFold"
      ;;
    a)
      echo "Ok, run settings within the folder with this experimentalID will be deleted, and the script processes as normal."
      rm -r "$simFoldDefaultFiles"
      rm -r "$simFoldOut"
      rm -r "$simFoldErr"
      ;;
    *)
      echo "You entered any string other than y or n, so the script has stopped anyway."
      exit
  esac
fi

## set and create additional basic folders ##
mkdir -p "$simFold"
mkdir -p "$simFoldDefaultFiles"
mkdir -p "$simFoldOut"
mkdir -p "$simFoldErr"

## validate user-defined model-options ##
validateUserConfigValues "$siteDataSet" "$validatesiteDataSet" "siteDataSet"
validateUserConfigValues "$customSiteList" "$validateCustomSiteList" "customSiteList"
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$queue" "$validateQueue" "queue"
validateUserConfigValues "$emailNotification" "$validateEmailNotification" "emailNotification"
validateUserConfigValues "$specialFluxnetOutput" "$validateSpecialFluxnetOutput" "specialFluxnetOutput"
validateUserConfigValues "$flag_calc_soil_uptake" "$validateFlags" "flag_calc_soil_uptake"
validateUserConfigValues "$flag_calc_soil_moisture" "$validateFlags" "flag_calc_soil_moisture"
validateUserConfigValues "$flag_base_soil_model" "$validateFlags" "flag_base_soil_model"
validateUserConfigValues "$flag_tile_veg_dynamics" "$validateFlags" "flag_tile_veg_dynamics"
validateUserConfigValues "$flag_start_bareground" "$validateFlags" "flag_start_bareground"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_optimal_Nleaf" "$validateFlags" "flag_optimal_Nleaf"
validateUserConfigValues "$flag_optimal_allocation" "$validateFlags" "flag_optimal_allocation"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$useSiteSpecificAWC" "$validateUseSiteSpecificAWC" "useSiteSpecificAWC"

## site data-set: FLUXNET or CRUNCEPv5 ##
if [ $siteDataSet == "f" ]  # if site is set to FLUXNET else CRUNCEPv5
 then
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET/"
 else 
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv5/"  # 64 sites, direct copy of 'sv_sites' but symlinks removed 
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/CRUNCEPv5/"
fi

## site list: default or custom ##
if [ $customSiteList == "yes" ]
 then
  # if $customSiteListFileName is NOT empty use the 
  if [ ! -z $customSiteListFileName ]
   then
    siteListFileName="$customSiteListFileName"
   else
    siteListFileName="my_site_list.dat"
  fi
  # test if file exists and has some content (no content-quality check!)
  if [ ! -e "$siteListFileName" ] # if file not exists
   then
    echo "EXIT! ""$siteListFileName"" does not exist! You set the option customSiteList to yes."
    exit
   else 
    if [ ! -s "$siteListFileName" ] # if file is empty
     then 
    echo "EXIT! ""$siteListFileName"" has no content, i.e. is empty! You set the option customSiteList to yes."
    exit
    fi
  fi
  # test if forcing data exist in the site-list-folder specified, for all sites specified in the custom site list
  siteList_lineNumber=$(wc -l "$siteListFileName" | awk '{print $1}')
  for (( a=2; a<=$siteList_lineNumber; a++ ))
    do
      siteName=$(sed -n "$a"p "$siteListFileName" | awk '{print $1}')  # -n  == quiet/silent
      if ! ls -la $siteDataFold | grep $siteName > /dev/null  # if forcing file not found
        then
          echo "EXIT! The site ""$siteName"" from your custom site_list.dat does not exist in the siteDataSet folder (f [FLUXNET], c [CRUNCEPv5]): ""$siteDataSet"" !"
          rm -r "$simFold"  # deleting the already created simulation folder
          exit
      fi
  done
  # copy custom site-list and set link to file accordingly
  cp "$siteListFileName" "$simFoldDefaultFiles"
  siteList="$simFoldDefaultFiles""$siteListFileName"
 else # the default case: use the fluxnet or cruncepv5 all_site_list.dat
  siteListFileName="all_site_list.dat"
  siteList="$siteDataFold""$siteListFileName"   # 2do --> may be customized to my_sitelist... later
fi

## validate experimental ID for alphanumeric only ##
if [ ! -z $(echo $jobName | tr -d "[:alnum:]" | tr -d "_" | tr -d "-") ]  # if "removing alphanum '_' '-' from jobName returns an non-empty string" then exit
 then
  echo "EXIT! Please use an experimental ID that is alphanumeric, i.e. without any special characters."
  exit
fi

## validate options linked with specialFluxnetOutput ##
if [ $specialFluxnetOutput == "yes" ]
 then
    if [ $siteDataSet != "f" ] 
      then
        echo "EXIT! You may use a FLUXNET site when you enable special fluxnet-simulation output."
        exit
    fi
fi

## validate allocated memory ##
lowestMem=10000       # 10 Mb
largestMem=50000000   # 50 Gb
# test memory empty or containing characters other then digits
if [ -z $memory ]  # returns TRUE if unset or empty
 then
  echo "EXIT! No memory defined. Please define the memory per job you think is useful."
  exit  
 else
  if [ ! -z $(echo $memory | tr -d "[:digit:]") ] # if "removing all digits" is not empty, i.e. other characters are included
   then
    echo "EXIT! The memory parameter contains characters other than digits. This is not working. Please use just numbers."
    exit
  fi
fi
# test against useful min and max values
if [ "$memory" -lt "$lowestMem" ]
 then
  echo "EXIT! Please increase memory, you allocated less than 10Mb; units is in Kb."
  exit
 else
  if [ "$memory" -gt "$largestMem" ]
   then
    echo "EXIT! Please decrease memory, you allocated more than 50Gb. Note, memory is allocated per job."
    exit
  fi
fi

## validate allocated number of CPUs ##
lowestCpu=5      # 5   parallel jobs
largestCpu=100   # 100 parallel jobs
# test memory empty or containing characters other then digits
if [ -z $numberCpuMax ]  # returns TRUE if unset or empty
 then
  echo "EXIT! No numberCpuMax values defined. Please define the number of CPUs that should run in parallel at the Cluster."
  exit  
 else
  if [ ! -z $(echo $numberCpuMax | tr -d "[:digit:]") ] # if "removing all digits" is not empty, i.e. other characters are included
   then
    echo "EXIT! The numberCpuMax parameter contains characters other than digits. This is not working. Please use just numbers."
    exit
  fi
fi
# test against useful min and max values
if [ "$numberCpuMax" -lt "$lowestCpu" ]
 then
  echo "Note: This is a Cluster. You may consider running more than 4 simulations in parallel."
 else
  if [ "$numberCpuMax" -gt "$largestCpu" ]
   then
    echo "EXIT! Please decrease the number of parallel jobs, you aim at using more than 100 CPUs in parallel."
    exit
  fi
fi

## convert $emailNotification into a bsup option ##
emailSwitch=""
if [ "$emailNotification" == "yes" ]
  then
    emailSwitch="-N"
fi

## test if namelist binary siteList R-Script exists ##
if [ ! -e "$slmBinary" ]     # exit if non-existing
  then
	echo "EXIT! binary is missing"
	exit
fi
if [ ! -e "$slmNamelist" ]   # exit if non-existing
  then
	echo "EXIT! namelist is missing"
	exit
fi
if [ ! -e "$siteList" ]      # exit if non-existing
  then
	echo "EXIT! site-list-data file is missing"
	exit
fi	
if [ ! -e "$slmRscript" ]    # exit if non-existing
  then
	echo "EXIT! the R Script is missing"
	exit
fi	

## validate $singleSiteName if user wants to run SLM only at one specific site ##
# validate if $singleSiteName is not empty
if [ ! -z "$singleSiteName" ]
 then
  allSites=$(cat "$siteList" | awk '{print $1}')
  validateUserConfigValues "$singleSiteName" "$allSites" "singleSiteName"
fi

## pass info about files to run_slm.sh via an extra file ("$simFoldDefaultFiles""config.txt") ##
echo "$siteList"              >  "$simFoldDefaultFiles""config.txt"
echo "$siteDataFold"          >> "$simFoldDefaultFiles""config.txt"
echo "$slmModel"              >> "$simFoldDefaultFiles""config.txt"
echo "$specialFluxnetOutput"  >> "$simFoldDefaultFiles""config.txt"
echo "$siteDataObsFold"       >> "$simFoldDefaultFiles""config.txt"
echo "$useSiteSpecificAWC"    >> "$simFoldDefaultFiles""config.txt"
#echo "$"    >> "$simFoldDefaultFiles""config.txt"
#echo "$"    >> "$simFoldDefaultFiles""config.txt"
#echo "$"    >> "$simFoldDefaultFiles""config.txt"

## copy namelist binary R-Script to simulation folder ##
cp "$slmBinary"     "$simFoldDefaultFiles"
cp "$slmNamelist"   "$simFoldDefaultFiles"
cp "$slmRscript"    "$simFoldDefaultFiles"
if [ $specialFluxnetOutput == "yes" ]; then 
  cp "$slmRscript4fluxfile" "$simFoldDefaultFiles"
  cp "$slmRscript4fluxdiag" "$simFoldDefaultFiles"
fi

## syncronize all data with file system ##
sync

## define a function that manipulates the namelist && change namelist according to options ##
# model: slm_model_name = land plant soil canopy
# 2 do: add other options like optimization

 # this function is copied from run_slm.sh without any modification
 # function: pass the variable_name ($1) & the new_value ($2)
 # $slmNamelist refers to the namelist file name, assuming the namelist file is located in the current folder !
function changeNamelistValue {
oldValue=$(grep -w $1 "$simFoldDefaultFiles""$slmNamelistFilename" | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 "$simFoldDefaultFiles""$slmNamelistFilename" | cut -d : -f 1)          # get line number of variable $1
sed -i "$lineNum"s/$oldValue/$2/ "$simFoldDefaultFiles""$slmNamelistFilename"                # replace the old value in line $lineNum with $2
sync
return 0                                                         # get return value with: return_value=$?
}
 # examples:
 #changeNamelistValue simulation_length_number 300
 #changeNamelistValue output_interval_flux monthly    # quotation marks not needed
 #changeNamelistValue flag_calc_soil_moisture .FALSE.
 #changeNamelistValue soil_theta_prescribe 122.7

changeNamelistValue slm_model_name "$slmModel"
if [ "$specialFluxnetOutput" == "yes" ]
  then
    changeNamelistValue fluxnet_output .TRUE.
  else
    changeNamelistValue fluxnet_output .FALSE.
fi
if [ $flag_calc_soil_uptake == T ]
  then
    changeNamelistValue flag_calc_soil_uptake ".TRUE."
  else
    changeNamelistValue flag_calc_soil_uptake ".FALSE."
fi
if [ $flag_calc_soil_moisture == T ]
  then
    changeNamelistValue flag_calc_soil_moisture ".TRUE."
  else
    changeNamelistValue flag_calc_soil_moisture ".FALSE."
fi
if [ $flag_base_soil_model == T ]
  then
    changeNamelistValue flag_base_soil_model ".TRUE."
  else
    changeNamelistValue flag_base_soil_model ".FALSE."
fi
if [ $flag_tile_veg_dynamics == T ]
  then
    changeNamelistValue flag_tile_veg_dynamics ".TRUE."
  else
    changeNamelistValue flag_tile_veg_dynamics ".FALSE."
fi
if [ $flag_start_bareground == T ]
  then
    changeNamelistValue flag_start_bareground ".TRUE."
  else
    changeNamelistValue flag_start_bareground ".FALSE."
fi
if [ $flag_t_jmax_acclimation == T ]
  then
    changeNamelistValue flag_t_jmax_acclimation ".TRUE."
  else
    changeNamelistValue flag_t_jmax_acclimation ".FALSE."
fi
if [ $flag_optimal_Nleaf == T ]
  then
    changeNamelistValue flag_optimal_Nleaf ".TRUE."
  else
    changeNamelistValue flag_optimal_Nleaf ".FALSE."
fi
if [ $flag_optimal_allocation == T ]
  then
    changeNamelistValue flag_optimal_allocation ".TRUE."
  else
    changeNamelistValue flag_optimal_allocation ".FALSE."
fi
if [ $flag_optimal_Nfraction == T ]
  then
    changeNamelistValue flag_optimal_Nfraction ".TRUE."
  else
    changeNamelistValue flag_optimal_Nfraction ".FALSE."
fi
changeNamelistValue simulation_length_number "$simulation_length_number"
changeNamelistValue output_interval_pool "$output_interval_pool"
changeNamelistValue output_interval_flux "$output_interval_flux"
changeNamelistValue output_start_first_day_year "$output_start_first_day_year"
changeNamelistValue output_end_last_day_year "$output_end_last_day_year"
changeNamelistValue soil_awc_prescribe "$soil_awc_prescribe"
changeNamelistValue soil_theta_prescribe "$soil_theta_prescribe"


## start job with simulations by starting the job(array) -> run_slm.sh ##
if [ -z $singleSiteName ]  # returns TRUE if unset or empty -> $siteList is used
 then # job array #
  array="yes"   # needed by cleanup_slm_hpc_run.sh
  # get number of lines in the file (number of sites stored in the $siteList = $siteNumber -1)
  siteNumber=$(wc -l "$siteList" | awk '{print $1}')
  # compile job array argument 
  jobArray=$(echo "$jobName"[2-"$siteNumber"]%"$numberCpuMax")
  # submit job array starting at 2 because first line is the header
  bsub -q "$queue" -M "$memory" -J "$jobArray" -o "$simFoldOut"out.%I.txt -e "$simFoldErr"err.%I.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles"
 else # job with individual site #
  array="no"   # needed by cleanup_slm_hpc_run.sh
  bsub -q "$queue" -M "$memory" -J "$jobName" -o "$simFoldOut"out.txt -e "$simFoldErr"err.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles" "$singleSiteName"
fi


## cleanup temporary files and folders via cleanup_slm_hpc_run.sh ##
# this job is pending until teh above started $jobName has ended
sleep 2s  # just some waiting time to avoid any interference with the above started job
cleanupJobName="$jobName""Cleanup"
bsub -w "ended($jobName)" -q "$queue" -J "$cleanupJobName" -o "$simFold"cleanup.eo "$emailSwitch" ./cleanup_slm_hpc_run.sh "$jobName" "$array" "$singleSiteName"

## end prep_run_slm.sh ##
