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

## folder name where data will be stored is: /Net/Groups/BGI/scratch/USER/slm/exp_'experimentalID'
# set something useful for you (allowed only: alphanumeric _ -); if empty ( experimentalID="" ) current date is used
experimentalID=""
## specify SLM model -> the parameter 'slm_model_name' in the namelist.slm is changed accordingly (so, no need to manually change it in the namelist)
# land (plant+soil), plant (necessary soil parts prescribed), soil (uses soil-forcing generated from plant), canopy
slmModel="land"


##>> sites <<###########################################################################################
## select set of sites: FLUXNET or CRUNCEPv7
# forcing data and all_site_list.dat for the two set of sites: /Net/Groups/BGI/data/Model/SLM1/input/point/ 
# options: fluxnet cruncep face
forcingDataSet="fluxnet"
# forcing with or without a transient (like spinup) phase before the actual simulation starts
# options: static transient
forcingMode="static"
# set the specific set of sites of the forcingDataSet you want to use
#  fluxnet:  (currently no options)
#  cruncep: all cue spp
#  face:     (currently no options)
siteSelection="cue"
# all sites or just some sites: specify a local custom site list if you want only some sites"
# T (True) or F (False) [large caps]
customSiteList="F"
# default siteList fileName: "my_site_list.dat" >> must be located in same folder as this script (usually "contrib/site_level/")
# you can specify a different file name other than my_site_list.dat; leave blank to use the default file name
customSiteListFileName=""
## run a single site: provide the Site-ID 
# note: this option must be empty (i.e., singleSiteName="" ) to run multiple sites
# Site-ID must be consistent with the forcingDataSet option above
# get the Site-ID from the 1st column of the particular "site list .dat"
# examples: [FLUX] "DE-Hai" "ZA-Kru" [CRUN] "BHS_01" "ILO_05" [FACE] "DK-AMB" "OR-ELE"
singleSiteName=""


##>> special runs <<####################################################################################
##(1) run Fluxnet simulations with data output at each timestep for the last complete forcing period at the end of the "simulation_length_number" you specify below
# T (True) or F (False) [large caps]
# when set to "T" simulations and extra scripts will be run by the run.slm.sh for each of the fluxnet sites you specify above (i.e., all sites, customSiteList, one site: singleSiteName)
# scripts are: convert_slm2FLUXNET.R & generate_FLUXNET_diagnostics.R
specialFluxnetOutput="F"

##(2) compare simulated model results of CRUNCEP sites with Global Forest Data Base data (/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_21092016.csv)
# T (True) or F (False) [large caps]
# when set to "T" options are set to: forcingDataSet="cruncep", output_interval_pool/flux="daily", paramSensitivitySim="F", flag_generate_soil_forcing="F", specialFluxnetOutput="F"
# the cleanup_slm_hpc_run.sh starts two R scripts: sim_cruncep_vs_gfdb-process_daily.R & sim_cruncep_vs_gfdb-create_figs.R
# you can run simulations for all/some CRUNCEPv7 sites of the "cue" set (cue_gfdb_sites*.dat) ; for only some sites (customSiteList="F/T"), or for just one site (singleSiteName="Site-ID")
compareSimCRUNCEPvsGFDB="F"
# years of spin-up before output is saved
spinupLengthGFDB=300

##(3) parameter sensitivity simulations
# you can run one site from your customSiteList (customSiteList is set to "T" automatically if paramSensitivitySim="T")
# T (True) or F (False) [large caps]
paramSensitivitySim="F"
# please provide the line number of the site you want to use in your "customSiteListFileName" >> line 2 to ...  (there is a header at line 1)
paramSensitivitySiteNumber=5 

##(4) generate forcing for stand alone soil model
# create additional output (veg%litterfall + some soilphysics) used as input of the (standalone) soil model
# T (True) or F (False) [large caps]
flag_generate_soil_forcing="F"
# years of spin-up before output is generated and saved
spinupLengthSoilForcing=200 
# number of years of forcing to be generated and saved after spin-up
yearsSoilForcing=1 


##>> namelist.slm options <<############################################################################
# flags: T (True) or F (False) [large caps]
# plant-soil model options
soil_moisture_scheme="layer"    # prescribed bucket layer
sb_model_scheme="simple_0d"     # simple_0d simple_1d jsm 
sb_nloss_scheme="dynamic"       # none fixed dynamic
flag_sb_prescribe_nh4="F"       # set nh4solute to a constant value
flag_sb_prescribe_no3="F"       # set no3solute to a constant value
flag_sb_prescribe_po4="T"       # set po4solute to a constant value
n_deposition_scheme="none"      # options: none constant dynamic (nhx_srf_down & noy_srf_down)
p_deposition_scheme="none"      # options: none constant dynamic (p_srf_down)

# BNF & Mycorrhiza
bnf_scheme="standard"           # none fixed standard optimal resistance
flag_mycorrhiza="F"             # TRUE: the other flag_myco* can be used with TRUE | FALSE: the other flag_myco* are set to FALSE by this script
flag_myco_inorganic="F"         # allows mycorrhiza to access the inorganic N pools
flag_myco_organic="F"           # allows mycorrhiza to access the organic/residual N pool
flag_myco_priming="F"           # present mycorrhiza enhance decomposition of residual pool (priming)

# canopy model options
flag_t_resp_acclimation="T"     # use the respiration temperature acclimation factor
flag_t_jmax_acclimation="T"     # acclimation of optimum temperature for Jmax
flag_optimal_Nfraction="F"      # optimal N fraction
ncanopy=10                      # number of canopy layers
canopy_layer_scheme="fapar"     # canopy layer depth calculation options: standard fapar

# plant growth model options
flag_optimal_Nleaf="F"          # optimal N leaf
flag_optimal_allocation="F"     # optimal allocation
# vegetation dynamics options
flag_tile_veg_dynamics="T"      # does also set "start_bareground"
flag_stand_harvest="F"          # if TRUE: l_do_stand_replacing_harvest is set to TRUE in SLM on the 1st timestep of the 1st day of the year "PlantYear" in the site_list file (or never, if that was before 1901)

# soil model options
soil_depth=5.7                  # meters

# simulation length & output
simulation_length_number=30     # default years for forcingMode: static = 30 | transient = 115 [transient: 01.01.1901 -> 31.12.2015 without spinup!]
output_interval_pool="weekly"   # timestep daily weekly monthly yearly
output_interval_flux="weekly"   # timestep daily weekly monthly yearly
output_start_first_day_year=1   # 1 to ... [must be lower than or equal with option: output_end_last_day_year]
output_end_last_day_year=30     # equal with or larger than option: output_start_first_day_year

# options for the forcingMode transient (only evaluated if forcingMode == transient)
transient_spinup_start_year=1901      # default 1901; start year of spinup period
transient_spinup_end_year=1930        # default 1930; last  year of spinup period
transient_spinup_years=500            # default 500; overall years of spinup repeating "spinup period" multiple times
transient_simulation_start_year=1901  # default 1901; simulation start with forcing at this year and ends after simulation_length_number

# forcing input
flag_oneYearForcing="F"         # T (True) or F (False): whether to use (repeatedly) only the 1st year of the climate forcing input file
soil_awc_prescribe=300.0        # float !
useSiteSpecificAWC="F"          # T (True) or F (False) [large caps] -> set "T" to use AWC as in site-list file (soil_awc_prescribe is used otherwise)
soil_theta_prescribe=1.0        # float !
flag_read_dC13="F"              # either use forcing dC13 or the prescribed default value
flag_read_dC14="F"              # either use forcing dC14 or the prescribed default value


##>> details about the HPC, if unsure stick with the defaults <<################################
queue="all"                  # available queues: all BTM io mpi_large  (for information on the queues run at pc026/28: bqueues -l)
emailNotification="no"       # set either to yes or no (small caps)
memory=1000                  # memory per job you may need [mega byte]; here provided in MB but recalculated to KB which is what OpenLava needs
numberCpuMax=80              # maximum number of simulations started in parallel


##>> flag_generate_soil_forcing options <<############################################################################
## these option-values apply specifically to the generation of soil forcing ##########################################
if [ $flag_generate_soil_forcing == "T" ]
 then
  echo " Simulation runs for soil_forcing creation:"
  echo " flag_generate_soil_forcing == T, applying special settings to some options"
  echo " please run the simulations (prep_run_slm.sh) for FLUXNET, and CRUNCEPv7 with siteSelection cue & spp"
  echo " "
  echo " after all jobs have finished, please run additional_scripts/copy_soil_forcing.sh"
  echo "   (this script updates the soil forcing for FLUXNET and CRUNCEPv7, "
  echo "    and renames the experimentalID but leaves the renamed simulation folders in your scratch)"
  echo " "
  experimentalID="create_soil_forc_""$forcingDataSet"  # the copy_soil_forcing.sh refers to this specific experimentalID
  if [ $forcingDataSet == "cruncep" ]
   then
    experimentalID="create_soil_forc_""$forcingDataSet""_""$siteSelection"
  fi
  forcingMode="static"
  slmModel="land"

  customSiteList="F"
  specialFluxnetOutput="F"
  compareSimCRUNCEPvsGFDB="F"
  paramSensitivitySim="F"

  soil_moisture_scheme="layer"
  sb_model_scheme="simple_0d"
  sb_nloss_scheme="dynamic"
  n_deposition_scheme="none"
  p_deposition_scheme="none"

  bnf_scheme="standard"
  flag_mycorrhiza="F"
  
  flag_t_resp_acclimation="T"
  flag_t_jmax_acclimation="T"
  flag_optimal_Nfraction="F"
  ncanopy=10 
  canopy_layer_scheme="fapar"
  flag_optimal_Nleaf="F"
  flag_optimal_allocation="F"
  flag_tile_veg_dynamics="T" 
  soil_depth=5.7  

  flag_read_dC13="F"
  flag_read_dC14="F"

  simulation_length_number=${spinupLengthSoilForcing}
  # forcing output is created per timestep independent of the below output_interval_pool/flux settings (i.e., write_soil_forcing() is called each timestep) !
  output_interval_pool="weekly"
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
while getopts ":e:f:c:m:s:r:y:o:b:z:n:h" opt; do
  case "$opt" in
  "e")
    echo "   experimentalID=""$OPTARG"
    experimentalID="$OPTARG"
    ;;
  "f")
    case "$OPTARG" in
    "fl")
      echo "   forcingDataSet=fluxnet"
      forcingDataSet="fluxnet"
      ;;
    "cr")
      echo "   forcingDataSet=cruncep"
      forcingDataSet="cruncep"
      ;;
    "fa")
      echo "   forcingDataSet=face"
      forcingDataSet="face"
      ;;
    esac
    ;;
  "c")
    echo "   customSiteList=""$OPTARG"
    customSiteList="$OPTARG"
    ;;
  "m")
    echo "   slmModel=""$OPTARG"
    slmModel="$OPTARG"
    ;;
  "s")
    case "$OPTARG" in
    "s0d")
      echo "   sb_model_scheme=simple_0d (only one soil layer)"
      sb_model_scheme="simple_0d"
      ;;
    "s1d")
      echo "   sb_model_scheme=simple_1d (multiple soil layers)"
      sb_model_scheme="simple_1d"
      ;;
    "jsm")
      echo "   sb_model_scheme=jsm (full-featured soil model)"
      sb_model_scheme="jsm"
      ;;
    esac
    ;;
  "y")
    echo "   simulation_length_number=""$OPTARG"
    echo "   output from 1 to ""$OPTARG"
    simulation_length_number="$OPTARG"
    output_start_first_day_year=1
    output_end_last_day_year="$OPTARG"
    ;;
  "r")
    case "$OPTARG" in
    "st")
      echo "   forcingMode=static"
      forcingMode="static"
      ;;
    "tr")
      echo "   forcingMode=transient"
      forcingMode="transient"
      ;;
    esac
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
  "b")
    case "$OPTARG" in
    "none")
      echo "   bnf_scheme=none"
      bnf_scheme="none"
      ;;
    "fix")
      echo "   bnf_scheme=fixed"
      bnf_scheme="fixed"
      ;;
    "std")
      echo "   bnf_scheme=standard"
      bnf_scheme="standard"
      ;;
    "opt")
      echo "   bnf_scheme=optimal"
      bnf_scheme="optimal"
      ;;
    "res")
      echo "   bnf_scheme=resistance"
      bnf_scheme="resistance"
      ;;
    esac
    ;;
  "z")
    case "$OPTARG" in
    "T")
      echo "   flag_mycorrhiza=""$OPTARG"
      echo "     mycorrhiza are active!"
      flag_mycorrhiza="$OPTARG"
      ;;
    "F")
      echo "   flag_mycorrhiza=""$OPTARG"
      echo "     Mycorrhiza are not active at all!"
      flag_mycorrhiza="$OPTARG"
      ;;
    esac
    ;;
  "n")
    echo "   your custom siteList filename is: ""$OPTARG"
    if [ $customSiteList == "F" ]
     then
      echo "     Anyway, you did not enable the use of this custom file (customSiteList is set to F), to enable you could use the option: -c T"
    fi
    customSiteListFileName="$OPTARG"
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
    echo " -f set the forcingDataSet"
    echo "    Arguments to -f:"
    echo "       fl: FLUXNET"
    echo "       cr: CRUNCEPv7"
    echo "       fa: FACE sites"
    echo " -c set the customSiteList"
    echo "    Arguments to -c:"
    echo "       T: using the custom site list"
    echo "       F: not using the custom site list"
    echo " -m set slmModel"
    echo "    Arguments to -m:"
    echo "       land plant soil canopy"
    echo " -s set the soil model"
    echo "    Arguments to -s:"
    echo "       s0d: simple_0d"
    echo "       s1d: simple_1d"
    echo "       jsm: jsm"
    echo " -y set simulation length"
    echo "    Arguments to -y:"
    echo "       any integer value"
    echo "       note: this sets simulation output accordingly (1-x)"
    echo " -r set forcingMode (namelist: forcing_mode)"
    echo "    Arguments to -r:"
    echo "       st: static"
    echo "       tr: transient"
    echo " -o set the optimality scheme"
    echo "    Arguments to -o:"
    echo "       no: no optimality"
    echo "       nf: using only optimal Nfraction"
    echo "       nl: using only optimal Nleaf"
    echo "       al: using only optimal allocation"
    echo " -b set the N uptake scheme"
    echo "    Arguments to -b:"
    echo "       none: no BNF"
    echo "       fix : fixed BNF"
    echo "       std : standard (includes costs)"
    echo "       opt : optimal"
    echo "       res : resistance"
    echo " -z set mycorrhiza on/off (flag_mycorrhiza)"
    echo "    Arguments to -z:"
    echo "       T: TRUE  - mycorrhiza are simulated"
    echo "       F: FALSE - mycorrhiza are NOT simulated"
    echo " -n define name of your custom site list"
    echo "    Any filename is possible, give the full name."
    echo "       e.g.: my_site_list.dat"
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
validateForcingDataSet="fluxnet cruncep face"
validateForcingMode="static transient"
validateSiteSelection="all cue spp"
validateSlmModel="land plant soil canopy"
validateBNFScheme="none fixed standard optimal resistance"
validateMycScheme="no am ecm"
validateNDepositionScheme="none constant dynamic"
validatePDepositionScheme="none constant dynamic"
validateCanopyLayerScheme="standard fapar"
validateSoilMoistureScheme="prescribed bucket layer"
validateSBNLossScheme="none fixed dynamic"
validateSBModelScheme="simple_0d simple_1d jsm"
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
slmLctlibFolder="../../data/"
slmLctlibFilename="lctlib_slm_nlct8.def"
slmParamSensiFolder="../parameter_list/"
slmParamListDefaultFilename="param.listD"
slmParamTableFilename="param_table.txt"
slmParamListDefault="$slmParamSensiFolder""param.listD"
slmParamTable="$slmParamSensiFolder""param_table.txt"
slmNamelist="$slmNamelistFolder""$slmNamelistFilename"            ## just one namelist for all models
slmLctlib="$slmLctlibFolder""$slmLctlibFilename"                  ## Landcovertype library for SLM with 8 PFTs
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
validateUserConfigValues "$forcingDataSet" "$validateForcingDataSet" "forcingDataSet"
validateUserConfigValues "$forcingMode" "$validateForcingMode" "forcingMode"
validateUserConfigValues "$siteSelection" "$validateSiteSelection" "siteSelection"
validateUserConfigValues "$customSiteList" "$validateFlags" "customSiteList"
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$queue" "$validateQueue" "queue"
validateUserConfigValues "$emailNotification" "$validateEmailNotification" "emailNotification"
validateUserConfigValues "$specialFluxnetOutput" "$validateFlags" "specialFluxnetOutput"
validateUserConfigValues "$bnf_scheme" "$validateBNFScheme" "bnf_scheme"
validateUserConfigValues "$flag_mycorrhiza" "$validateFlags" "flag_mycorrhiza"
validateUserConfigValues "$flag_myco_inorganic" "$validateFlags" "flag_myco_inorganic"
validateUserConfigValues "$flag_myco_organic" "$validateFlags" "flag_myco_organic"
validateUserConfigValues "$flag_myco_priming" "$validateFlags" "flag_myco_priming"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$soil_moisture_scheme" "$validateSoilMoistureScheme" "soil_moisture_scheme"
validateUserConfigValues "$sb_model_scheme" "$validateSBModelScheme" "sd_model_scheme"
validateUserConfigValues "$sb_nloss_scheme" "$validateSBNLossScheme" "sb_nloss_scheme"
validateUserConfigValues "$flag_sb_prescribe_nh4" "$validateFlags" "flag_sb_prescribe_nh4"
validateUserConfigValues "$flag_sb_prescribe_no3" "$validateFlags" "flag_sb_prescribe_no3"
validateUserConfigValues "$flag_sb_prescribe_po4" "$validateFlags" "flag_sb_prescribe_po4"
validateUserConfigValues "$flag_tile_veg_dynamics" "$validateFlags" "flag_tile_veg_dynamics"
validateUserConfigValues "$flag_t_resp_acclimation" "$validateFlags" "flag_t_resp_acclimation"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_optimal_Nleaf" "$validateFlags" "flag_optimal_Nleaf"
validateUserConfigValues "$flag_optimal_allocation" "$validateFlags" "flag_optimal_allocation"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$flag_generate_soil_forcing" "$validateFlags" "flag_generate_soil_forcing"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$flag_read_dC13" "$validateFlags" "flag_read_dC13"
validateUserConfigValues "$flag_read_dC14" "$validateFlags" "flag_read_dC14"
validateUserConfigValues "$useSiteSpecificAWC" "$validateFlags" "useSiteSpecificAWC"
validateUserConfigValues "$paramSensitivitySim" "$validateFlags" "paramSensitivitySim"
validateUserConfigValues "$compareSimCRUNCEPvsGFDB" "$validateFlags" "compareSimCRUNCEPvsGFDB"


## apply setting for comparison simulated model results of CRUNCEP sites with Global Forest Data Base ##
# this may overwrite some options specifies above !
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then
  gfdbFile="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_cue_sites_01082017.csv"
  forcingDataSet="cruncep"
  forcingMode="static"
  siteSelection="cue"
  paramSensitivitySim="F"
  flag_generate_soil_forcing="F"
  specialFluxnetOutput="F"
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
    echo " "
    echo " You started simulations to compare CRUNCEP with GFDB data using a spin-up time of ""$spinupLengthGFDB"" years."
    echo " The comparison will be based on daily SLM output of "$(($output_end_last_day_year - $output_start_first_day_year + 1))" years after the spin-up."
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


## site data-set: FLUXNET / CRUNCEPv7 / FACE ##
case "$forcingDataSet" in
"fluxnet")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/"
  #siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET2/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET/"
  ;;
"cruncep")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv7/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/CRUNCEPv7/"
  ;;
"face")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FACE/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FACE/"
  ;;
esac


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
  echo ""
  echo " You aim at parameter sensitivity simulations for one specific site from your customSiteList."
  echo " Simulations will be done for the parameter combinations you specified in contrib/parameter_list/param_table.txt"
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
      # for CRUNCEPv7 sites one forcing file is used for multiple sites differentiated with a number e.g. SITE_01, remove number for site test
      if [ $forcingDataSet == "cruncep" ]
       then
        siteName=$(cut -c1-3  <<< "$siteName")
      fi
      # test if a forcing file for this site is available
      if ! ls -la $siteDataFoldForcing | grep $siteName > /dev/null  # if forcing file not found
        then
          echo "EXIT! The site ""$siteName"" from your custom site_list.dat does not exist in the forcingDataSet folder of: ""$forcingDataSet"" !"
          rm -r "$simFold"  # deleting the already created simulation folder
          exit
      fi
  done
  # copy custom site-list and set link to file accordingly
  cp "$siteListFileName" "$simFoldDefaultFiles"
  siteList="$simFoldDefaultFiles""$siteListFileName"
else # the default case: $customSiteList == "F"
  ## select the correct site information file (a.k.a. "all_site_list.dat")
  ## file names differ across forcingDataSets & forcingModes
  case "$forcingDataSet" in
  "fluxnet")
    siteListFileName="all_sites_list.dat"
    ;;
  "cruncep")
    case "$siteSelection" in
    "all")
      if [ $forcingMode == "static" ]
       then
        siteListFileName="all_gfdb_sites_1984-2013_list.dat"
      else
        siteListFileName="all_gfdb_sites_1901-2015_list.dat"
      fi
      ;;
    "cue")
      if [ $forcingMode == "static" ]
       then
        siteListFileName="cue_gfdb_sites_1984-2013_list.dat"
      else
        siteListFileName="cue_gfdb_sites_1901-2015_list.dat"
      fi
      ;;
    "spp")
      if [ $forcingMode == "static" ]
       then
        siteListFileName="spp_sites_1984-2013_list.dat"
      else
        siteListFileName="spp_sites_1901-2015_list.dat"
      fi
      ;;
    esac
    ;;
  "face")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="all_sites_1996-2008_list.dat"
    else
      siteListFileName="all_sites_1901-2015_list.dat"
    fi
    ;;
  esac
  siteList="$siteDataFold""$siteListFileName"
fi


## evaluate options for the forcingMode == transient
if [ $forcingMode == "transient" ]
 then
   echo " "
   echo " You set the forcingMode to transient, with a spinup time of ""$transient_spinup_years"" years that uses forcing input from year ""$transient_spinup_start_year"" to ""$transient_spinup_end_year""."
   echo " The subsequent simulation may cover ""$simulation_length_number"" years, starting with forcing data from year ""$transient_simulation_start_year""."

   if [ $transient_spinup_start_year != "1901" ]
    then
     echo " You set the transient_spinup_start_year to ""$transient_spinup_start_year"", which is currently not possible. Sorry."
     echo " Only the default value (1901) is working. Setting transient_spinup_start_year=1901."
     transient_spinup_start_year=1901
   fi
   if [ $transient_simulation_start_year != "1901" ]
    then
     echo " You set the transient_simulation_start_year to ""$transient_simulation_start_year"", which is currently not possible. Sorry."
     echo " Only the default value (1901) is working. Setting transient_simulation_start_year=1901."
     transient_simulation_start_year=1901
   fi
   # test calc length of forcing period used for spinup
   transient_test_spinup_forcing_period=$((transient_spinup_end_year - transient_spinup_start_year + 1))
   if [ $transient_test_spinup_forcing_period -lt 6 ] || [ $transient_test_spinup_forcing_period -gt 49 ]
    then
     echo " You set the length of forcing period repeatedly used for spinup to ""$transient_test_spinup_forcing_period""."
     echo " This is not recommended. Please use a period length (transient_spinup_end_year - transient_spinup_start_year) between 5 and 50 years."
     transient_spinup_end_year=$((transient_spinup_start_year + 29))
     echo " The new (automagically corrected) value of transient_spinup_end_year is ""$transient_spinup_end_year"" setting the period length to the default of 30 years."
   fi
   if [ $transient_spinup_years -lt 10 ]
    then
     echo " You set the spinup length to ""$transient_spinup_years"", which quite low and not very useful."
     echo " Assuming a simple typo, the spinup length is set to 100 years (the default would be 500 years though)."
     transient_spinup_years=100
   fi
   if [ $forcingDataSet == "fluxnet" ]
    then
     echo " Please note: the transient forcingMode works best with CRUNCEPv7 and FACE sites, but might produce questionable results with FLUXNET2 sites."
     echo "  You selected FLUXNET2 sites. This is just a warning, none of your options will be changed."
   fi
 else  
 # $forcingMode == "static"
 # set transient* parameter to default values
  transient_spinup_start_year=1901
  transient_spinup_end_year=1930
  transient_spinup_years=500
  transient_simulation_start_year=1901
fi


## an additional identifier (ID) for static and transient forcing files -> passed to (via config file) and used in run_slm.sh
forcingModeFileID="0"
case "$forcingDataSet" in
"fluxnet")
  forcingModeFileID="_"
  ;;
"cruncep")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_1984"
  else
    forcingModeFileID="_1901"
  fi
  ;;
"face")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_199"
  else
    forcingModeFileID="_190"
  fi
  ;;
esac


## validate experimental ID for alphanumeric only ##
if [ ! -z $(echo $jobName | tr -d "[:alnum:]" | tr -d "_" | tr -d "-") ]  # if removing alphanum _ - from jobName returns an non-empty string" then exit
 then
  echo "EXIT! Please use an experimental ID that is alphanumeric, i.e. without any special characters (_ and - are ok though)."
  exit
fi


## set mycorrhiza flag_myco_* flags to FALSE if mycorrhiza are disabled overall
if [ $flag_mycorrhiza == "F" ]
 then
  flag_myco_inorganic="F"
  flag_myco_organic="F"
  flag_myco_priming="F"
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
      else
        echo " "
        echo "You just started simulations to generate FLUXNET diagnostics using the last full forcing period after ""$simulation_length_number"" years."
        #if $singleSiteName is not empty
        if [ ! -z "$singleSiteName" ]
         then
           echo "Simulation is done for one site: ""$singleSiteName"
         else
           echo "Simulations are done for all sites listed in: ""$siteListFileName"
        fi
    fi
    # check site data set
    if [ $forcingDataSet != "fluxnet" ] 
      then
        echo "EXIT! You might want to use (a) FLUXNET site(s) when enabling special FLUXNET-Simulation output :) ."
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
if [ ! -e "$slmLctlib" ]   # exit if non-existing
  then
	echo "EXIT! lctlib is missing"
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
echo "$forcingDataSet"                              >> "$simFoldDefaultFiles""config.txt"
echo "$forcingModeFileID"                           >> "$simFoldDefaultFiles""config.txt"
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

## copy namelist binary R-Scripts to simulation folder (i.e., at the scratch drive) ##
cp "$slmBinary"     "$simFoldDefaultFiles"
cp "$slmNamelist"   "$simFoldDefaultFiles"
cp "$slmLctlib"     "$simFoldDefaultFiles"
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

# change values in the namelist copied to the simulation folder at the scratch drive
changeNamelistValue slm_model_name "$slmModel"
changeNamelistValue forcing_mode "$forcingMode"
changeNamelistValue bnf_scheme "$bnf_scheme" 
changeNamelistValue n_deposition_scheme "$n_deposition_scheme"
changeNamelistValue p_deposition_scheme "$p_deposition_scheme"
changeNamelistValue ncanopy "$ncanopy"
changeNamelistValue canopy_layer_scheme "$canopy_layer_scheme" 
changeNamelistValue soil_moisture_scheme "$soil_moisture_scheme" 
changeNamelistValue sb_model_scheme "$sb_model_scheme" 
changeNamelistValue sb_nloss_scheme "$sb_nloss_scheme"
if [ $flag_mycorrhiza == T ]
  then
    changeNamelistValue flag_mycorrhiza ".TRUE."
  else
    changeNamelistValue flag_mycorrhiza ".FALSE."
fi
if [ $flag_myco_inorganic == T ]
  then
    changeNamelistValue flag_myco_inorganic ".TRUE."
  else
    changeNamelistValue flag_myco_inorganic ".FALSE."
fi
if [ $flag_myco_organic == T ]
  then
    changeNamelistValue flag_myco_organic ".TRUE."
  else
    changeNamelistValue flag_myco_organic ".FALSE."
fi
if [ $flag_myco_priming == T ]
  then
    changeNamelistValue flag_myco_priming ".TRUE."
  else
    changeNamelistValue flag_myco_priming ".FALSE."
fi
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
if [ $flag_tile_veg_dynamics == T ]
  then
    changeNamelistValue flag_tile_veg_dynamics ".TRUE."
  else
    changeNamelistValue flag_tile_veg_dynamics ".FALSE."
fi
if [ $flag_stand_harvest == T ]
  then
    changeNamelistValue flag_stand_harvest ".TRUE."
  else
    changeNamelistValue flag_stand_harvest ".FALSE."
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

if [ $flag_read_dC13 == T ]
  then
    changeNamelistValue flag_read_dC13 ".TRUE."
  else
    changeNamelistValue flag_read_dC13 ".FALSE."
fi
if [ $flag_read_dC14 == T ]
  then
    changeNamelistValue flag_read_dC14 ".TRUE."
  else
    changeNamelistValue flag_read_dC14 ".FALSE."
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
changeNamelistValue transient_spinup_start_year "$transient_spinup_start_year"
changeNamelistValue transient_spinup_end_year "$transient_spinup_end_year"
changeNamelistValue transient_spinup_years "$transient_spinup_years"
changeNamelistValue transient_simulation_start_year "$transient_simulation_start_year"
  

## start job with simulations by starting the job(array) -> run_slm.sh ##
if [ -z $singleSiteName ]  # returns TRUE if unset or empty -> $siteList is used
 then # job array #
  array="T"   # needed by cleanup_slm_hpc_run.sh
  # get number of lines in the file (number of sites stored in the $siteList = $siteNumber -1 and starts at line 2, see jobArray below)
  siteNumber=$(wc -l "$siteList" | awk '{print $1}')
  # overwrite siteNumber with number of parameter combinations for parameter sensitivity
  if [ $paramSensitivitySim == "T" ]; then
    siteNumber=$(($paramSensi_NumberCombinations +1)) # 2 to number+1
  fi
  # compile job array argument 
  jobArray=$(echo "$jobName"[2-"$siteNumber"]%"$numberCpuMax")
  # submit job array starting at 2 because first line is the header
  bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobArray" -o "$simFoldOut"out.%I.txt -e "$simFoldErr"err.%I.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles"
else # job with individual site #
  array="F"   # needed by cleanup_slm_hpc_run.sh
  bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobName" -o "$simFoldOut"out.txt -e "$simFoldErr"err.txt "$emailSwitch" ./run_slm.sh "$jobName" "$simDefaultFiles" "$singleSiteName"
fi


## cleanup temporary files and folders via cleanup_slm_hpc_run.sh ##
# this job is pending until the above started $jobName has ended
sleep 1s  # just some waiting time to avoid any interference with the above started job
cleanupJobName="$jobName""Cleanup"
cleanQueue="BTM all"
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then cleanQueue="io BTM"; fi
bsub -w "ended($jobName)" -q "$cleanQueue" -J "$cleanupJobName" -o "$simFold"cleanup.eo "$emailSwitch" ./cleanup_slm_hpc_run.sh "$jobName" "$array" "$singleSiteName"



## end prep_run_slm.sh ##
