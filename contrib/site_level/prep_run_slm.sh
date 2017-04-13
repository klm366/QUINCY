#!/bin/bash

## 
##
## > prepare and start run_slm.sh script to run SLM simulations across multiple sites
## > the prep_run_slm.sh & run_slm-sh should be located in trunk/contrib/site_level/ (same for any branch)
## > you can run this script from any folder, i.e. you have not to cd into contrib/site_level/ 
##   (this is because this script knows about its location, due to the command: 'cd "$( dirname "$0" )"')
## 

########################################################################################################
####### simulation specific options you might want to change ###########################################
########################################################################################################

# folder name where data will be stored is: /Net/Groups/BGI/scratch/USER/slm/exp_'experimentalID'
experimentalID=""      # set something useful for you (allowed only: alphanumeric _ -); if empty ( experimentalID="" ) current date is used
# specify SLM model -> the parameter 'slm_model_name' in the namelist.slm is changed accordingly (so, no need to manually change it in the namelist)
slmModel="land"             # land, plant, soil are currently available


##>> sites <<###########################################################################################
# select sites: FLUXNET or CRUNCEPv5
siteDataSet="f"           # f: FLUXNET or c: CRUNCEPv5  (use small caps!)
# want a custom site list? -> default name is: "my_site_list.dat"; must be located in same folder as this script (usually "contrib/site_level/")
customSiteList="F"        # T (True) or F (False) [large caps] -> the sites in your file must be constistent with the option siteDataSet f/c
customSiteListFileName="my_fluxnet_site_list.dat"    # you can specify a different file name other than my_site_list.dat; leave blank to use the default file name
# run only a single site? then provide a site name, else write: singleSiteName=""
singleSiteName=""         # one site name, e.g. "DE-Bay"; must be consistent with the parameter value: siteDataSet f/c; overwrites customSiteList


##>> special runs <<####################################################################################
# (1) run Fluxnet simulations with data output at each timestep for the last complete forcing period (see land_model_driver.f90)
specialFluxnetOutput="F"    # T (True) or F (False) [large caps]

# (2) compare simulated model results of CRUNCEP sites with Global Forest Data Base data (/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_21092016.csv)
# started by the cleanup_slm_hpc_run.sh using two R scripts: sim_cruncep_vs_gfdb-process_daily.R & sim_cruncep_vs_gfdb-create_figs.R
compareSimCRUNCEPvsGFDB="F"  # T (True) or F (False) [large caps] -> when "T" defaults are: siteDataSet="c", output_interval_pool/flux="daily", paramSensitivitySim="F"
spinupLengthGFDB=300         # years of spin-up before output is saved

# (3) parameter sensitivity simulations
# you can run one site from your customSiteList (customSiteList is set to "T" automatically if paramSensitivitySim="T")
# please provide the line number of the site you want to use in your "customSiteListFileName"
paramSensitivitySim="F"         # T (True) or F (False) [large caps]
paramSensitivitySiteNumber=5    # line 2 to ...  (there is a header at line 1)

# (4) generate forcing for stand alone soil model
flag_generate_soil_forcing="F"  # create additional output (veg%litterfall + some soilphysics) used as input of the (standalone) soil model
spinupLengthSoilForcing=200       # total length of run to produce soil forcing
yearsSoilForcing=1              # number of years of forcing to be generated    

##>> namelist.slm options <<############################################################################
# flags: T (True) or F (False) [large caps]
# plant-soil model options
n_uptake_scheme="standard"      # prescribed standard std_plus_costs optimal resistance
p_uptake_scheme="prescribed"    # prescribed standard 
soil_moisture_scheme="bucket"   # prescribed bucket layer
sb_model_scheme="one_pool"      # one_pool four_pool jsm 
flag_sb_prescribe_nh4="F"       # set nh4solute to a constant value
flag_sb_prescribe_no3="F"       # set no3solute to a constant value
flag_sb_prescribe_po4="F"       # set po4solute to a constant value
flag_n_deposition="F"           # enable/disable T/F the use of nhx_deposition & noy_deposition
flag_calc_soil_energy="F"
# canopy model options
flag_t_resp_acclimation="T"
flag_t_jmax_acclimation="T"
flag_optimal_Nfraction="F"
flag_heterogeneous_canopy="T"
ncanopy=20                          # number of canopy layers
canopy_layer_scheme="standard"      # canopy layer depth calculation options: standard fapar
# plant growth model options
flag_optimal_Nleaf="F"
flag_optimal_allocation="F"
# vegetation dynamics options
flag_tile_veg_dynamics="T"          # does also set "start_bareground"
# soil model options
soil_depth=5.7                  # meters

# simulation length & output
simulation_length_number=30     # 1 to ...500 ? [simulation_length_unit defaults to "y" == years]
output_interval_pool="weekly"   # timestep daily weekly monthly yearly
output_interval_flux="weekly"   # timestep daily weekly monthly yearly
output_start_first_day_year=1   # 1 to ... [lower than or equal with output_end_last_day_year ]
output_end_last_day_year=30     # equal with or larger than output_start_first_day_year

# forcing input
flag_oneYearForcing="F"         # T (True) or F (False): whether to use (repeatedly) only the 1st year of the climate forcing input file
soil_awc_prescribe=300.0        # float !
useSiteSpecificAWC="F"          # T (True) or F (False) [large caps] -> set "T" to use AWC as in site-list file (soil_awc_prescribe is used otherwise)
soil_theta_prescribe=1.0      # float !


##>> details about the HPC, if unsure stick with these defaults <<################################
queue="BTM"                  # available queues: all io BTM mpi_large
fastNodes="F"                # if T only the nodes iof iog ioh ioi ioj are used (works only with queue BTM and io) | T (True) or F (False) [large caps]
emailNotification="no"       # set either to yes or no (small caps)
memory=1000                  # memory per job you may need [mega byte]; here provided in MB but recalculated to KB which is what OpenLava needs
numberCpuMax=80              # maximum number of simulations started in parallel


##>> flag_generate_soil_forcing options <<############################################################################
## these option-values apply specifically to the generation of soil forcing ##########################################
if [ $flag_generate_soil_forcing == "T" ]
 then
  echo " flag_generate_soil_forcing == T, applying special settings to some options"
  echo " please run the simulations (prep_run_slm.sh) for both FLUXNET and CRUNCEPv5"
  echo " afterwards, please run additional_scripts/copy_soil_forcing.sh when all jobs have finished"
  echo " "
  experimentalID="create_soil_forc_""$siteDataSet"  # the copy_soil_forcing.sh refers to this specific experimentalID
  slmModel="plant"

  customSiteList="F"
  specialFluxnetOutput="F"
  compareSimCRUNCEPvsGFDB="F"
  paramSensitivitySim="F"

  n_uptake_scheme="standard"
  p_uptake_scheme="standard"
  soil_moisture_scheme="layer"
  sb_model_scheme="one_pool"
  flag_calc_soil_energy="T"

  flag_t_resp_acclimation="T"
  flag_t_jmax_acclimation="T"
  flag_optimal_Nfraction="F"
  flag_heterogeneous_canopy="T"
  ncanopy=20 
  canopy_layer_scheme="standard"
  flag_optimal_Nleaf="F"
  flag_optimal_allocation="F"
  flag_tile_veg_dynamics="T" 
  soil_depth=5.7  

  simulation_length_number=${spinupLengthSoilForcing}
  output_interval_pool="weekly"  # write_soil_forcing() is called per timestep independent of the output_interval_.. settings
  output_interval_flux="weekly"
  let fd=spinupLengthSoilForcing-yearsSoilForcing+1
  output_start_first_day_year=${fd}
  output_end_last_day_year=${spinupLengthSoilForcing}

  flag_oneYearForcing="F"
  soil_awc_prescribe=300.0
  useSiteSpecificAWC="F"
  soil_theta_prescribe=1.0
fi


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"


## deal with options and arguments passed to the script ##
while getopts ":e:s:c:m:y:o:n:h" opt; do
  case "$opt" in
  "e")
    echo "   experimentalID=""$OPTARG"
    experimentalID="$OPTARG"
    ;;
  "s")
    if [ "$OPTARG" == "f" ]
     then
      echo "   siteDataSet=f (FLUXNET)"
     else
      echo "   siteDataSet=c (CRUNCEPv5)"
    fi
    siteDataSet="$OPTARG"
    ;;
  "c")
    echo "   customSiteList=""$OPTARG"
    customSiteList="$OPTARG"
    ;;
  "m")
    echo "   slmModel=""$OPTARG"
    slmModel="$OPTARG"
    ;;
  "y")
    echo "   simulation_length_number=""$OPTARG"
    echo "   output from 1 to ""$OPTARG"
    simulation_length_number="$OPTARG"
    output_start_first_day_year=1
    output_end_last_day_year="$OPTARG"
    ;;
  "o")
    case "$OPTARG" in
    "no")
      echo "   no optimization"
      flag_optimal_Nfraction="F"
      flag_optimal_Nleaf="F"
      flag_optimal_allocation="F"
      ;;
    "nf")
      echo "   using only optimal Nfraction"
      flag_optimal_Nfraction="T"
      flag_optimal_Nleaf="F"
      flag_optimal_allocation="F"
      ;;
    "nl")
      echo "   using only optimal Nleaf"
      flag_optimal_Nfraction="F"
      flag_optimal_Nleaf="T"
      flag_optimal_allocation="F"
      ;;
    "al")
      echo "   using only optimal allocation"
      flag_optimal_Nfraction="F"
      flag_optimal_Nleaf="F"
      flag_optimal_allocation="T"
      ;;
    esac
    ;;
  "n")
    case "$OPTARG" in
    "pre")
      echo "   n_uptake_scheme=prescribed"
      n_uptake_scheme="prescribed"
      ;;
    "std")
      echo "   n_uptake_scheme=standard"
      n_uptake_scheme="standard"
      ;;
    "stdp")
      echo "   n_uptake_scheme=std_plus_costs"
      n_uptake_scheme="std_plus_costs"
      ;;
    "opt")
      echo "   n_uptake_scheme=optimal"
      n_uptake_scheme="optimal"
      ;;
    "res")
      echo "   n_uptake_scheme=resistance"
      n_uptake_scheme="resistance"
      ;;
    esac
    ;;
  \?)
    echo "     Invalid option: -$OPTARG"
    echo "     Use -h for help."
    ;;
  :)
    echo "Option -$OPTARG requires an argument."
    echo "EXIT !"
    exit 1
    ;;
  "h")  # help #
    echo "Valid options are:"
    echo " -e set the experimentalID"
    echo " -s set the siteDataSet"
    echo "    Arguments to -s:"
    echo "       f: FLUXNET"
    echo "       c: CRUNCEPv5"
    echo " -c set the customSiteList"
    echo "    Arguments to -c:"
    echo "       T: using the custom site list"
    echo "       F: not using the custom site list"
    echo " -m set slmModel"
    echo "    Arguments to -m:"
    echo "       land plant soil"
    echo " -y set simulation length"
    echo "    Arguments to -y:"
    echo "       any integer value"
    echo "       note: this sets simulation output accordingly (1-x)"
    echo " -o set the optimality scheme"
    echo "    Arguments to -o:"
    echo "       no: no optimality"
    echo "       nf: using only optimal Nfraction"
    echo "       nl: using only optimal Nleaf"
    echo "       al: using only optimal allocation"
    echo " -n set the N uptake scheme"
    echo "    Arguments to -n:"
    echo "       pre : prescribed"
    echo "       std : standard"
    echo "       stdp: standard plus transformation costs"
    echo "       opt : optimal"
    echo "       res : resistance"
    exit 1
    ;;
  esac
done


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
validateSlmModel="land plant soil"  # canopy
validateNUptakeScheme="prescribed standard std_plus_costs optimal resistance"
validatePUptakeScheme="prescribed standard"
validateCanopyLayerScheme="standard fapar"
validateSoilMoistureScheme="prescribed bucket layer"
validateSBModelScheme="one_pool four_pool jsm"
validateQueue="all io BTM mpi_large"
validateEmailNotification="yes no"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validateOptionForcing="constant generate readfile"

## set and create basic folders and files ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
mkdir -p "$bgiScratchUserFold"
slmBinary="../../x86_64-intel/bin/land.x"                  ## just one binary for all models
slmNamelistFolder="../namelist/"
slmNamelistFilename="namelist.slm"
slmParamSensiFolder="../parameter_list/"
slmParamListDefaultFilename="param.listD"
slmParamTableFilename="param_table.txt"
slmParamListDefault="$slmParamSensiFolder""param.listD"
slmParamTable="$slmParamSensiFolder""param_table.txt"
slmNamelist="$slmNamelistFolder""$slmNamelistFilename"            ## just one namelist for all models
slmRscript="../r_script/postprocessing_""$slmModel"".R"           ## just the R-Script for postprocessing output data
slmRscript4fluxfile="../r_script/convert_slm2FLUXNET.R"           ## just the R-Script for converting output for FLUXNET 
slmRscript4fluxdiag="../r_script/generate_FLUXNET_diagnostics.R"  ## just the R-Script for FLUXNET diagnostics 
slmRscriptCompCruncepGfdb1="../r_script/sim_cruncep_vs_gfdb-process_daily.R"  ## just the R-Script for postprocessing simulation data for CRUNCEP vs GFDB
slmRscriptCompCruncepGfdb2="../r_script/sim_cruncep_vs_gfdb-create_figs.R"    ## just the R-Script for creating comparision figures CRUNCEP vs GFDB

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
simParamSensiFiles="parameter_sensitivity_files/"
simFoldDefaultFiles="$simFold""$simDefaultFiles"
simFoldParamSensiFiles="$simFold""$simParamSensiFiles"
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


## validate user-defined model-options ##
validateUserConfigValues "$siteDataSet" "$validatesiteDataSet" "siteDataSet"
validateUserConfigValues "$customSiteList" "$validateFlags" "customSiteList"
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$queue" "$validateQueue" "queue"
validateUserConfigValues "$emailNotification" "$validateEmailNotification" "emailNotification"
validateUserConfigValues "$specialFluxnetOutput" "$validateFlags" "specialFluxnetOutput"
validateUserConfigValues "$n_uptake_scheme" "$validateNUptakeScheme" "n_uptake_scheme"
validateUserConfigValues "$p_uptake_scheme" "$validatePUptakeScheme" "p_uptake_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$soil_moisture_scheme" "$validateSoilMoistureScheme" "soil_moisture_scheme"
validateUserConfigValues "$sb_model_scheme" "$validateSBModelScheme" "sd_model_scheme"
validateUserConfigValues "$flag_sb_prescribe_nh4" "$validateFlags" "flag_sb_prescribe_nh4"
validateUserConfigValues "$flag_sb_prescribe_no3" "$validateFlags" "flag_sb_prescribe_no3"
validateUserConfigValues "$flag_sb_prescribe_po4" "$validateFlags" "flag_sb_prescribe_po4"
validateUserConfigValues "$flag_n_deposition" "$validateFlags" "flag_n_deposition"
validateUserConfigValues "$flag_calc_soil_energy" "$validateFlags" "flag_calc_soil_energy"
validateUserConfigValues "$flag_tile_veg_dynamics" "$validateFlags" "flag_tile_veg_dynamics"
validateUserConfigValues "$flag_t_resp_acclimation" "$validateFlags" "flag_t_resp_acclimation"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_optimal_Nleaf" "$validateFlags" "flag_optimal_Nleaf"
validateUserConfigValues "$flag_optimal_allocation" "$validateFlags" "flag_optimal_allocation"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$flag_heterogeneous_canopy" "$validateFlags" "flag_heterogeneous_canopy"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$flag_generate_soil_forcing" "$validateFlags" "flag_generate_soil_forcing"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$useSiteSpecificAWC" "$validateFlags" "useSiteSpecificAWC"
validateUserConfigValues "$paramSensitivitySim" "$validateFlags" "paramSensitivitySim"
validateUserConfigValues "$compareSimCRUNCEPvsGFDB" "$validateFlags" "compareSimCRUNCEPvsGFDB"
validateUserConfigValues "$fastNodes" "$validateFlags" "fastNodes"


## apply setting for comparison simulated model results of CRUNCEP sites with Global Forest Data Base ##
# this may overwrite some options specifies above !
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then
  gfdbFile="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_21092016.csv"
  siteDataSet="c"
  paramSensitivitySim="F"
  let simulation_length_number=spinupLengthGFDB+simulation_length_number
  output_interval_pool="daily"
  output_interval_flux="daily"
  let output_start_first_day_year=spinupLengthGFDB+output_start_first_day_year
  let output_end_last_day_year=spinupLengthGFDB+output_end_last_day_year     
  if [ ! -e "$gfdbFile" ] # if file not exists
   then
    echo "EXIT! The GFDB file does not exist! You set the option compareSimCRUNCEPvsGFDB to T."
    echo " File location & name should be: ""$gfdbFile"
    exit
   else
    echo " Note: compareSimCRUNCEPvsGFDB was set to T."
    echo "  You aim for daily SLM simulation output across "$(($output_end_last_day_year - $output_start_first_day_year + 1))" years."
  fi
fi

## set and create additional basic folders (within the basic simulation folder) ##
mkdir -p "$simFold"
mkdir -p "$simFoldDefaultFiles"
mkdir -p "$simFoldOut"
mkdir -p "$simFoldErr"
if [ $paramSensitivitySim == "T" ]; then
  mkdir -p "$simFoldParamSensiFiles"
fi

## site data-set: FLUXNET or CRUNCEPv5 ##
if [ $siteDataSet == "f" ]  # if site is set to FLUXNET else CRUNCEPv5
 then
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET/"
 else 
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv5/"  # 64 sites, direct copy of 'sv_sites' but symlinks removed 
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/CRUNCEPv5/"
fi


## SLM model forcing: land/plant/canopy/soil 
case "$slmModel" in
"land")
  siteDataFoldForcing="$siteDataFold""met_forcing/"
  ;;
"plant")
  siteDataFoldForcing="$siteDataFold""met_forcing/"
  ;;
"canopy")
  siteDataFoldForcing="$siteDataFold""met_forcing/"
  ;;
"soil")
  siteDataFoldForcing="$siteDataFold""soil_forcing/"
  ;;
esac


## test site list: default or custom ##
# parameter sensitivity: set customSiteList to "T" and delete any name in singleSiteName (see the IF[] for starting the job very below)
if [ $paramSensitivitySim == "T" ]; then
  customSiteList="T"
  singleSiteName=""
fi
# test default and custom site list
if [ $customSiteList == "T" ]
 then
  # if option $customSiteListFileName is NOT empty use the 
  if [ ! -z $customSiteListFileName ]
   then
    siteListFileName="$customSiteListFileName"
   else
    siteListFileName="my_site_list.dat"
  fi
  # test if custom-sitelist-file exists and has some content (no content-quality check!)
  if [ ! -e "$siteListFileName" ] # if file not exists
   then
    echo "EXIT! ""$siteListFileName"" does not exist! You set the option customSiteList to T."
    exit
   else 
    if [ ! -s "$siteListFileName" ] # if file is empty
     then 
    echo "EXIT! ""$siteListFileName"" has no content, i.e. is empty! You set the option customSiteList to T."
    exit
    fi
  fi
  # test if forcing data exist in the site-list-folder specified, for all sites specified in the custom site list
  # test custom-sitelist only
  siteList_lineNumber=$(wc -l "$siteListFileName" | awk '{print $1}') # number of sites in the custom-sitelist
  for (( a=2; a<=$siteList_lineNumber; a++ )) # line number 2 to ... excluding the header at line 1
    do
      siteName=$(sed -n "$a"p "$siteListFileName" | awk '{print $1}')  # -n  == quiet/silent
      if ! ls -la $siteDataFoldForcing | grep $siteName > /dev/null  # if forcing file not found
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
  siteList="$siteDataFold""$siteListFileName"
fi

## validate experimental ID for alphanumeric only ##
if [ ! -z $(echo $jobName | tr -d "[:alnum:]" | tr -d "_" | tr -d "-") ]  # if removing alphanum _ - from jobName returns an non-empty string" then exit
 then
  echo "EXIT! Please use an experimental ID that is alphanumeric, i.e. without any special characters (_ and - are ok though)."
  exit
fi

## validate options linked with specialFluxnetOutput ##
if [ $specialFluxnetOutput == "T" ]
 then
    # check SLM model
    if [ $slmModel == "soil" ]
      then
        echo "EXIT! specialFluxnetOutput currently not working with SLM ""$slmModel"" model."
        rm -r "$simFold"  # deleting the already created simulation folder
        exit
    fi
    # check site data set
    if [ $siteDataSet != "f" ] 
      then
        echo "EXIT! You may use a FLUXNET site when you enable special FLUXNET-Simulation output :) ."
        exit
    fi
fi

## validate allocated memory ##
lowestMem=10       # 10 Mb
largestMem=50000   # 50 Gb
# test memory empty or containing characters other then digits
if [ -z $memory ]  # returns TRUE if unset or empty
 then
  echo "EXIT! No memory defined. Please define the memory per job you think is useful (and within a reasonable range)."
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
  echo "EXIT! Please increase memory, you allocated less than 10Mb. Note, memory is defined with unit Kilobyte [Kb]."
  exit
 else
  if [ "$memory" -gt "$largestMem" ]
   then
    echo "EXIT! Please decrease memory, you allocated more than 50Gb. Note, memory is allocated per job."
    exit
  fi
fi


## compile the -R argument (Reserved); specifically the "reserved memory"
reservedArg="rusage[mem="$memory"]"


## calculate memory (-M argument) from MB -> KB as needed by OpenLava
memory=$((memory * 1000))


## validate allocated number of CPUs ##
lowestCpu=5      # 5   parallel jobs
largestCpu=100   # 100 parallel jobs
# test numberCpuMax empty or containing characters other then digits
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
  echo "Note: This is a High Performance Cluster. You may consider running more than 4 simulations in parallel."
 else
  if [ "$numberCpuMax" -gt "$largestCpu" ]
   then
    echo "EXIT! Please decrease the number of parallel jobs, you aim at using more than 100 CPUs in parallel."
    echo "(Of course, technically you could use more than 100 CPUs. This limit is simply set to be nice with your colleagues.)"
    exit
  fi
fi

## convert $emailNotification into a bsup option ##
emailSwitch=""
if [ "$emailNotification" == "yes" ]
  then
    emailSwitch="-N"
fi

## set fastNodes option for BTM queue ##
if [ $fastNodes == "T" ]
  then
    if [ $queue == "BTM" ] || [ $queue == "io" ]  # BTM or io
      then
        fN="iof iog ioh ioi ioj"
      else
        echo "Note: the fast nodes [iof iog ioh ioi ioj] are available only for the BTM and io queue. I switch to fastNodes = F."
        fastNodes="F"
    fi
fi

## test if namelist binary siteList R-Script (and param.listD + param_table.txt) exists ##
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
# in case of parameter sensitivity
if [ $paramSensitivitySim == "T" ]
  then
    # param.listD
    if [ ! -e "$slmParamListDefault" ]    # exit if non-existing
      then
    	echo "EXIT! the default parameter List (param.listD) is missing though needed for Parameter Sensitivity Simulations."
      echo " You may run the contrib/parameter_list/generate_parameter_table.sh to create this file."
    	exit
    fi
    # param_table.txt
    if [ ! -e "$slmParamTable" ]    # exit if non-existing
      then
    	echo "EXIT! the parameter Table (param_table.txt) is missing though needed for Parameter Sensitivity Simulations."
      echo " You may run the contrib/parameter_list/generate_parameter_table.sh to create the table."
    	exit
    fi	
fi

## validate $singleSiteName if user wants to run SLM only at one specific site ##
# validate if $singleSiteName is not empty
if [ ! -z "$singleSiteName" ]
 then
  allSites=$(cat "$siteList" | awk '{print $1}')
  validateUserConfigValues "$singleSiteName" "$allSites" "singleSiteName"
fi

## parameter sensitivity: get the number of different parameter combinations ##
if [ $paramSensitivitySim == "T" ]
  then
    paramSensi_NumberCombinations=$(wc -l "$slmParamTable" | awk '{print $1}')
    paramSensi_NumberCombinations=$(($paramSensi_NumberCombinations -1))  # not counting the header of the table
    # validate number of parameter combinations
    if [ $paramSensi_NumberCombinations -gt 994 ]
      then
        echo "EXIT !  You want to do ""$paramSensi_NumberCombinations"" simulations. That is too much. A number lower than 995 is necessary."
        exit
      else
        echo "The number of combinations of parameter values (i.e. simulations) is: ""$paramSensi_NumberCombinations"
    fi
fi

## pass info about files to run_slm.sh via an extra file ("$simFoldDefaultFiles""config.txt") ##
echo "$siteList"                                    >  "$simFoldDefaultFiles""config.txt"
echo "$siteDataFoldForcing"                         >> "$simFoldDefaultFiles""config.txt"
echo "$slmModel"                                    >> "$simFoldDefaultFiles""config.txt"
echo "$specialFluxnetOutput"                        >> "$simFoldDefaultFiles""config.txt"
echo "$siteDataObsFold"                             >> "$simFoldDefaultFiles""config.txt"
echo "$useSiteSpecificAWC"                          >> "$simFoldDefaultFiles""config.txt"
echo "$paramSensitivitySim"" ""$simParamSensiFiles" >> "$simFoldDefaultFiles""config.txt"
echo "$compareSimCRUNCEPvsGFDB"                     >> "$simFoldDefaultFiles""config.txt"
#echo "$"                          >> "$simFoldDefaultFiles""config.txt"
#echo "$"                          >> "$simFoldDefaultFiles""config.txt"

# pass extra info about parameter sensitivity simulations to folder
if [ $paramSensitivitySim == "T" ]; then
  echo "$paramSensitivitySiteNumber"                  >  "$simFoldParamSensiFiles""param_config.txt"
  echo "$slmParamListDefaultFilename"                 >>  "$simFoldParamSensiFiles""param_config.txt"
  echo "$slmParamTableFilename"                       >>  "$simFoldParamSensiFiles""param_config.txt"
  #echo "$"                  >>  "$simFoldParamSensiFiles""param_config.txt"
  #echo "$"                  >>  "$simFoldParamSensiFiles""param_config.txt"
fi

## copy namelist binary R-Scripts to simulation folder ##
cp "$slmBinary"     "$simFoldDefaultFiles"
cp "$slmNamelist"   "$simFoldDefaultFiles"
cp "$slmRscript"    "$simFoldDefaultFiles"
if [ $specialFluxnetOutput == "T" ]; then 
  cp "$slmRscript4fluxfile" "$simFoldDefaultFiles"
  cp "$slmRscript4fluxdiag" "$simFoldDefaultFiles"
fi
if [ $paramSensitivitySim == "T" ]; then
  cp "$slmParamListDefault" "$simFoldParamSensiFiles"
  cp "$slmParamTable"       "$simFoldParamSensiFiles"
fi
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then
  cp "$slmRscriptCompCruncepGfdb1" "$simFoldDefaultFiles"
  cp "$slmRscriptCompCruncepGfdb2" "$simFoldDefaultFiles"
fi

## syncronize all data with the file system ##
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
 #changeNamelistValue soil_theta_prescribe 1.

changeNamelistValue slm_model_name "$slmModel"
changeNamelistValue n_uptake_scheme "$n_uptake_scheme"
changeNamelistValue p_uptake_scheme "$p_uptake_scheme"
changeNamelistValue ncanopy "$ncanopy"
changeNamelistValue canopy_layer_scheme "$canopy_layer_scheme" 
changeNamelistValue soil_moisture_scheme "$soil_moisture_scheme" 
changeNamelistValue sb_model_scheme "$sb_model_scheme" 
if [ $flag_sb_prescribe_nh4 == T ]
  then
    changeNamelistValue flag_sb_prescribe_nh4 ".TRUE."
  else
    changeNamelistValue flag_sb_prescribe_nh4 ".FALSE."
fi
if [ $flag_sb_prescribe_no3 == T ]
  then
    changeNamelistValue flag_sb_prescribe_no3 ".TRUE."
  else
    changeNamelistValue flag_sb_prescribe_no3 ".FALSE."
fi
if [ $flag_sb_prescribe_po4 == T ]
  then
    changeNamelistValue flag_sb_prescribe_po4 ".TRUE."
  else
    changeNamelistValue flag_sb_prescribe_po4 ".FALSE."
fi
if [ $flag_n_deposition == T ]
  then
    changeNamelistValue flag_n_deposition ".TRUE."
  else
    changeNamelistValue flag_n_deposition ".FALSE."
fi
if [ "$specialFluxnetOutput" == "T" ]
  then
    changeNamelistValue fluxnet_output .TRUE.
  else
    changeNamelistValue fluxnet_output .FALSE.
fi
if [ $paramSensitivitySim == T ]
  then
    changeNamelistValue set_parameter_values_from_file ".TRUE."
  else
    changeNamelistValue set_parameter_values_from_file ".FALSE."
fi
if [ $flag_calc_soil_energy == T ]
  then
    changeNamelistValue flag_calc_soil_energy ".TRUE."
  else
    changeNamelistValue flag_calc_soil_energy ".FALSE."
fi
if [ $flag_tile_veg_dynamics == T ]
  then
    changeNamelistValue flag_tile_veg_dynamics ".TRUE."
  else
    changeNamelistValue flag_tile_veg_dynamics ".FALSE."
fi
if [ $flag_t_resp_acclimation == T ]
  then
    changeNamelistValue flag_t_resp_acclimation ".TRUE."
  else
    changeNamelistValue flag_t_resp_acclimation ".FALSE."
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
if [ $flag_heterogeneous_canopy == T ]
  then
    changeNamelistValue flag_heterogeneous_canopy ".TRUE."
  else
    changeNamelistValue flag_heterogeneous_canopy ".FALSE."
fi
if [ $flag_oneYearForcing == T ]
  then
    changeNamelistValue flag_climate_forcing_one_year ".TRUE."
  else
    changeNamelistValue flag_climate_forcing_one_year ".FALSE."
fi
if [ $flag_generate_soil_forcing == T ]
  then
    changeNamelistValue flag_generate_soil_forcing ".TRUE."
  else
    changeNamelistValue flag_generate_soil_forcing ".FALSE."
fi
changeNamelistValue simulation_length_number "$simulation_length_number"
changeNamelistValue output_interval_pool "$output_interval_pool"
changeNamelistValue output_interval_flux "$output_interval_flux"
changeNamelistValue output_start_first_day_year "$output_start_first_day_year"
changeNamelistValue output_end_last_day_year "$output_end_last_day_year"
changeNamelistValue soil_awc_prescribe "$soil_awc_prescribe"
changeNamelistValue soil_theta_prescribe "$soil_theta_prescribe"
changeNamelistValue soil_depth "$soil_depth"



## start job with simulations by starting the job(array) -> run_slm.sh ##
if [ -z $singleSiteName ]  # returns TRUE if unset or empty -> $siteList is used
 then # job array #
  array="T"   # needed by cleanup_slm_hpc_run.sh
  # get number of lines in the file (number of sites stored in the $siteList = $siteNumber -1)
  siteNumber=$(wc -l "$siteList" | awk '{print $1}')
  # overwrite siteNumber with number of parameter combinations for parameter sensitivity
  if [ $paramSensitivitySim == "T" ]; then
    siteNumber=$(($paramSensi_NumberCombinations +1)) # 2 to number+1
  fi
  # compile job array argument 
  jobArray=$(echo "$jobName"[2-"$siteNumber"]%"$numberCpuMax")
  # submit job array starting at 2 because first line is the header
  if [ $fastNodes == "T" ]
    then
      bsub -q "$queue" -m "$fN" -M "$memory" -R "$reservedArg" -J "$jobArray" -o "$simFoldOut"out.%I.txt -e "$simFoldErr"err.%I.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles"
    else
      bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobArray" -o "$simFoldOut"out.%I.txt -e "$simFoldErr"err.%I.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles"
  fi
 else # job with individual site #
  array="F"   # needed by cleanup_slm_hpc_run.sh
  if [ $fastNodes == "T" ]
    then
      bsub -q "$queue" -m "$fN" -M "$memory" -R "$reservedArg" -J "$jobName" -o "$simFoldOut"out.txt -e "$simFoldErr"err.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles" "$singleSiteName"
    else
      bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobName" -o "$simFoldOut"out.txt -e "$simFoldErr"err.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles" "$singleSiteName"
  fi
fi


## cleanup temporary files and folders via cleanup_slm_hpc_run.sh ##
# this job is pending until the above started $jobName has ended
sleep 2s  # just some waiting time to avoid any interference with the above started job
cleanupJobName="$jobName""Cleanup"
bsub -w "ended($jobName)" -q "$queue" -J "$cleanupJobName" -o "$simFold"cleanup.eo "$emailSwitch" ./cleanup_slm_hpc_run.sh "$jobName" "$array" "$singleSiteName"

## end prep_run_slm.sh ##
