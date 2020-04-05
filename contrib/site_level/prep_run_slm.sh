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
# land (plant+soil), plant (necessary soil parts prescribed), soil (uses soil-forcing generated from land), canopy
slmModel="land"


##>> sites <<###########################################################################################
## select the site-set:
# forcing data and site info *list.dat for all site-sets: /Net/Groups/BGI/data/Model/SLM1/input/point/ 
# options: fluxnet fluxnet3 cruncep face manip fert craine photacc spp1685 confin
forcingDataSet="fluxnet"
# forcing with or without a transient (like spinup) phase before the actual simulation starts
# options: static transient
forcingMode="static"
# set the specific set of sites of the forcingDataSet you want to use
#  cruncep: all cue spp
crunSites="cue"
# soil-forcing revision number
# revision of the trunk the soil-forcing was created with (to see whats avilable, check /Net/Groups/BGI/data/Model/SLM1/input/point/"forcingDataSet"/ )
# for static and transient soil-forcing (transient only for: fluxnet/3)
# provide only numbers, e.g., 2174  (not "r2174")
soilForcStatRevNum=""
soilForcTransRevNum="2174"
# all sites or just some sites: specify a local custom site list if you want only some sites"
# T (True) or F (False) [large caps]
customSiteList="F"
# default siteList fileName: "my_site_list.dat" >> must be located in same folder as this script (usually "contrib/site_level/")
# you can specify a different file name other than my_site_list.dat; leave blank to use the default file name
customSiteListFileName=""
## run a single site: provide the Site-ID 
# this option must be empty (i.e., singleSiteName="" ) to run multiple sites !
# Site-ID must be consistent with the forcingDataSet option above, i.e., fluxnet, fluxnet3, cruncep, face, manip, ...
# you can get the Site-ID from the 1st column of the particular site-list info file in /Net/Groups/BGI/data/Model/SLM1/input/point/"forcingDataSet"/
# examples: [FLUX/3] "DE-Hai" "FI-Hyy" "ZA-Kru"  [CRUN] "BHS_01" "HAI_01" "ILO_05"  
#           [FACE] "DK-AMB" "OR-ELE"  [MANIP] "ES-MaG" "ES-MaT"  [FERT] "HAR-CON_4" "HAR-HIG_4" "HAR-LOW_4" "HAR-CON_5"
#           [CRAINE] "CRA_01" "CRA_02" "CRA_18"  [PHOTACC] "ARC_01" "MON_02" "STR_01"  [SPP_1685] "BBR" "LUE" "VES"
#           [CONFIN] "FI-Hyy" "FI-Ken"
singleSiteName=""


##>> special runs <<####################################################################################
##(1) fluxnet/3 diagnostics output 
# run Fluxnet/3 simulations creating timestep output for the last full forcing period within the simulation runtime 'simulation_length_number'
# in addition output is created at the interval 'output_interval_pool' (must not be "timestep") from 'output_start_...' to 'output_end_...'
# T (True) or F (False) [large caps]
# when set to "T" simulations and extra scripts will be run by the run.slm.sh for each of the fluxnet/3 sites you specify above (i.e., all sites, customSiteList, one site: singleSiteName)
# scripts are: convert_slm2FLUXNET.R & generate_FLUXNET_diagnostics.R
specialFluxnetOutput="F"

##(2) cruncep vs gfdb comparison
# compare simulated model results of CRUNCEP sites with Global Forest Data Base data (/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_21092016.csv)
# T (True) or F (False) [large caps]
# when set to "T" options are set to: paramSensitivitySim="F", specialFluxnetOutput="F"
# output is created from 'output_start_...' to 'output_end_...' (typically 1986 to 2015), but see the output_interval_*_spinup options
# the output period may be set to output_interval_pool/flux="daily" so that the automatic model-data comparison works
# the cleanup_slm_hpc_run.sh starts the R scripts: sim_cruncep_vs_gfdb-process_daily.R & sim_cruncep_vs_gfdb-create_figs.R & sim_cruncep_all_vs_gfdb-create_figs.R
# you can run simulations for the site-sets 'all' and 'cue' CRUNCEPv7, or for only some/one site(s) of it (use customSiteList="F/T"; or singleSiteName="Site-ID")
compareSimCRUNCEPvsGFDB="F"
# years of spin-up (static and transient forcing) before output is saved (all output is used for simulation vs GFDB comparision)
# forcingMode static:    simulation_length_number is increased and 'output_start_...' & 'output_end_...' are shifted by the years defined as spinup
# forcingMode transient: transient_spinup_years=spinupLengthGFDB
spinupLengthGFDB=500

##(3) parameter sensitivity simulations
# you can run a set of sensitivity simulations for one site; T (True) or F (False) [large caps]
paramSensitivitySim="F"
# you can give absolute values ("F") or proportional changes within range [0 < x < 2] ("T")
paramChangesProportional="F"
# please provide the Site-ID   (see examples in the three lines above the option "singleSiteName" ~ line 61)
paramSensitivitySiteName="DE-Hai"


##>> default usecases <<################################################################################
## set the according 'usecase number' here (e.g. usecase="20")
## or run: prep_run_slm.sh -u 'usecase number' (e.g. prep_run_slm.sh -u 20 )
## see contrib/usecases/ for available usecases
usecase=""


##>> namelist.slm options <<############################################################################
# flags: T (True) or F (False) [large caps]
# plant-soil model options
sb_model_scheme="simple_0d"     # simple_0d simple_1d jsm 
sb_nloss_scheme="dynamic"       # none fixed dynamic
flag_sb_prescribe_nh4="F"       # set nh4solute to a constant value
flag_sb_prescribe_no3="F"       # set no3solute to a constant value
flag_sb_prescribe_po4="T"       # set po4solute to a constant value

# BNF & Mycorrhiza
bnf_scheme="standard"           # none fixed standard optimal resistance unlimited
flag_mycorrhiza="F"             # TRUE: the other flag_myco* can be used with TRUE | FALSE: the other flag_myco* are set to FALSE by this script
flag_mycorrhiza_org="F"         # allows mycorrhiza to access the organic/residual N pool
flag_mycorrhiza_prim="F"        # present mycorrhiza enhance decomposition of residual pool (priming)

# canopy model options
flag_t_resp_acclimation="T"     # use the respiration temperature acclimation factor
flag_t_jmax_acclimation="T"     # acclimation of optimum temperature for Jmax
flag_optimal_Nfraction="F"      # optimal N fraction
ncanopy=10                      # number of canopy layers
canopy_layer_scheme="fapar"     # canopy layer depth calculation options: standard fapar

# plant growth model options
biomass_alloc_scheme="fixed"    # fixed dynamic optimal
leaf_stoichom_scheme="fixed"    # fixed dynamic optimal
# vegetation dynamics options
veg_dynamics_scheme="population"  # none population cohort
cohort_harvest_interval=80        # number of years after which the cohort is harvested and re-planted, is used with veg_dynamics_scheme="cohort"
flag_stand_harvest="F"            # if TRUE: l_do_stand_replacing_harvest is set to TRUE in SLM on the 1st timestep of the 1st day of the year "PlantYear" in the site_list file (or never, if that was before 1901)

# soil model options
soil_depth=5.7                  # meter (5.7 to 9.5)
nsoil_layer=5                   # number of soil layers (default for SLM land: 5 | for SLM soil only available forcing: 15 | generally: min 5 & max 20)

# JSM
flag_sb_jsm_transport="F"         # ..
flag_sb_jsm_litter_input="F"      # ..
flag_sb_jsm_OM_sorption="F"       # ..

# ECA and P cycle
sb_adsorp_scheme="eca_full"       # eca_full; eca_part (only uptake, not for adsorption); eca_none

# simulation length & output
simulation_length_number=30       # default years for forcingMode: static = 30 | transient = 115 [transient: 01.01.1901 -> 31.12.2015 without spinup!]
output_interval_pool="weekly"     # timestep daily weekly monthly yearly
output_interval_flux="weekly"     # timestep daily weekly monthly yearly
output_start_first_day_year=1     # 1 to ... [must be lower than or equal with option: output_end_last_day_year]
output_end_last_day_year=30       # equal with or larger than option: output_start_first_day_year

# options for the forcingMode transient (only evaluated if forcingMode == transient)
transient_spinup_start_year=1901      # default 1901; start year of spinup period
transient_spinup_end_year=1930        # default 1930; last  year of spinup period
transient_spinup_years=500            # default 500; overall years of spinup repeating "spinup period" multiple times
transient_simulation_start_year=1901  # default 1901; simulation start with forcing at this year and ends after simulation_length_number
output_interval_pool_spinup="none"    # none daily weekly monthly yearly (only for transient forcingMode; if you'd need timestep please talk to Jan or Sönke)
output_interval_flux_spinup="none"    # none daily weekly monthly yearly (only for transient forcingMode; if you'd need timestep please talk to Jan or Sönke)
flag_transientSim_CO2const="F"        # T (True) or F (False): if TRUE the first value of CO2 (dC13 DC14) from 1901 is used as constant for the full simulation after spinup (i.e. no variation in CO2 dC13 DC14)

# forcing input
flag_oneYearForcing="F"         # T (True) or F (False): whether to use (repeatedly) only the 1st year of the climate forcing input file
soil_awc_prescribe=300.0        # float !
useSiteSpecificAWC="F"          # T (True) or F (False) [large caps] -> set "T" to use AWC as in site-list file (soil_awc_prescribe is used otherwise)
soil_theta_prescribe=1.0        # float !
flag_read_dC13="F"              # either use forcing dC13 or the prescribed default value
flag_read_DC14="F"              # either use forcing DC14 or the prescribed default value
n_deposition_scheme="none"      # options: none constant dynamic (nhx_srf_down & noy_srf_down)
p_deposition_scheme="none"      # options: none constant dynamic (p_srf_down)


##>> details about the HPC, if unsure stick with the defaults <<################################
queue="all"                  # available queues: all BTM io mpi_large io_perform (for information on the queues run at pc026/27/28: bqueues -l; io_perform selects only the new INTEL3 nodes: lshosts)
emailNotification="no"       # set either to yes or no (small caps)
memory=5000                  # memory per job you may need [mega byte]; here provided in MB (default: 1000) but recalculated to KB which is what OpenLava needs
numberCpuMax=160             # maximum number of simulations started in parallel (default is 160)




########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"


## deal with options and arguments passed to the script ##
# getopts
while getopts ":e:f:c:m:s:r:y:b:z:n:u:h" opt; do
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
    "f3")
      echo "   forcingDataSet=fluxnet3 (Berkley 2015)"
      forcingDataSet="fluxnet3"
      ;;
    "cr")
      echo "   forcingDataSet=cruncep"
      forcingDataSet="cruncep"
      ;;
    "fa")
      echo "   forcingDataSet=face"
      forcingDataSet="face"
      ;;
    "ma")
      echo "   forcingDataSet=manip"
      forcingDataSet="manip"
      ;;
    "fe")
      echo "   forcingDataSet=fert"
      forcingDataSet="fert"
      ;;
    "ca")
      echo "   forcingDataSet=craine"
      forcingDataSet="craine"
      ;;
    "ph")
      echo "   forcingDataSet=photacc"
      forcingDataSet="photacc"
      ;;
    "sp")
      echo "   forcingDataSet=spp1685"
      forcingDataSet="spp1685"
      ;;
    "co")
      echo "   forcingDataSet=confin"
      forcingDataSet="confin"
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
    "unl")
      echo "   bnf_scheme=unlimited"
      bnf_scheme="unlimited"
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
  "u")
    usecase="$OPTARG"
    echo "   usecase=""$OPTARG"
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
    echo "       f3: FLUXNET3 (Berkley 2015)"
    echo "       cr: CRUNCEPv7"
    echo "       fa: FACE sites"
    echo "       ma: MANIP tree and grassland sites"
    echo "       fe: FERT sites (fertilization)"
    echo "       ca: CRAINE sites (selected sites from doi:10.5061/dryad.v2k2607)"
    echo "       ph: PHOTACC sites (phot acclimation)"
    echo "       sp: SPP_1685 sites (for jsm simulations)"
    echo "       co: CONFIN (four finnish sites with conifers)"
    echo " -c set the customSiteList"
    echo "    Arguments to -c:"
    echo "       T: using the custom site list"
    echo "       F: not using the custom site list"
    echo " -n define name of your custom site list"
    echo "    Any filename is possible, give the full name."
    echo "       e.g.: my_site_list.dat"
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
    echo " -b set the N uptake scheme"
    echo "    Arguments to -b:"
    echo "       none: no BNF"
    echo "       fix : fixed BNF"
    echo "       std : standard (includes costs)"
    echo "       opt : optimal"
    echo "       res : resistance"
    echo "       unl : unlimited"
    echo " -z set mycorrhiza on/off (flag_mycorrhiza)"
    echo "    Arguments to -z:"
    echo "       T: TRUE  - mycorrhiza are simulated"
    echo "       F: FALSE - mycorrhiza are NOT simulated"
    echo " -u provide the two-digit number of the usecase; for available usecases see contrib/usecases/ "
    echo "    You can also create new custom usecase files (see the wiki for explanation)."
    exit 1
    ;;
  esac
done


## include usecase
if [ ! -z "$usecase" ]            # if $usecase is not empty
 then
   usecaseFileIncl=$( echo "../usecases/""$usecase""_"* )
   if [ -e "$usecaseFileIncl" ]   # if $usecase file is present
    then
      . "$usecaseFileIncl"        # include text
      echo " "
      echo " Applying settings for usecase: ""$usecase""  (usecase file: ""$usecaseFileIncl"")"
      echo " "
      # test whether the experimentalID is defined in the usecase file, 
      # if so:  append the revision number of the local working copy
      #         also used for soil-forcing creation (usecases 01 & 08)
      validateExpIDPresent="$(grep -w "experimentalID" "$usecaseFileIncl")"
      if [ ! -z "${validateExpIDPresent// }" ]; then  # -z test is explicitly for: the length of string is zero, here not zero
        # check if SVN information is available (exit status == 0), if so add the current revision number to the experimentalID
        svn info ./ >/dev/null 2>&1
        exitSVN=$?
        if [ "$exitSVN" == 0 ]; then # if exit status of svn info is zero, i.e., all good
          workingCopyRev=$(svn info ../../ | grep "Rev:" | awk '{print $4}')
          experimentalID="$experimentalID""_r""$workingCopyRev"
          echo " The script could retrieve the revision number of your local working copy."
          echo " The experimentalID is set to: ""$experimentalID"
          echo " "
        else  # if rev number could not be retrieved, append "0000"
          experimentalID="$experimentalID""_r0000"
        fi
      fi
   else 
     echo " "
     echo " You want to apply the settings of usecase number: ""$usecase"
     echo " This failed because the include file: ""$usecaseFileIncl"" is missing."
     echo " EXIT!"
     exit 1
  fi
else
  usecase="none"
fi


## detect and deal with the fast / veryfast / extrafast argument that might be passed to the script ##
# by using one of these arguments, essentially, the simulations are:
# i)  bound to the 'all' queue
# ii) the reserved memory per job is increased in a way that only a limited number of jobs (simulations) per node are possible
fastOpt="FALSE"
vfastOpt="FALSE"
efastOpt="FALSE"
fastArg="$(grep -w "fast" <<< "$@")"
vfastArg="$(grep -w "veryfast" <<< "$@")"
efastArg="$(grep -w "extrafast" <<< "$@")"
# if the list of arguments to the prep_run_slm.sh script does contain 'fast', i.e. $fastArgs is not empty
if [ ! -z "${fastArg// }" ]  # -z test is explicitly for: the length of string is zero
 then 
  fastOpt="TRUE"
  queue="all"
  memory=18000    # three jobs per all queue node | potentially 90 SLM jobs can run in parallel
  echo "  Info: you opted for the fast mode, simulations are done at the all queue with a max of 3 SLM jobs per node."
  echo " "
fi
# if the list of arguments to the prep_run_slm.sh script does contain 'veryfast', i.e. $vfastArgs is not empty
if [ ! -z "${vfastArg// }" ]  # -z test is explicitly for: the length of string is zero
 then 
  vfastOpt="TRUE"
  queue="all"
  memory=28000     # two jobs per all queue node | potentially 60 SLM jobs can run in parallel
  echo "  Info: you opted for the veryfast mode, simulations are done at the all queue with a max of 2 SLM jobs per node."
  echo " "
fi
# if the list of arguments to the prep_run_slm.sh script does contain 'extrafast', i.e. $efastArgs is not empty
if [ ! -z "${efastArg// }" ]  # -z test is explicitly for: the length of string is zero
 then 
  efastOpt="TRUE"
  queue="all"
  memory=58000     # one job per all queue node | note: the all queue consists of "only" 30 nodes
  echo "  Info: you opted for the extrafast mode, simulations are done at the all queue with only 1 SLM job per node."
  echo " "
fi


## prepare Validation of user-defined model-options ##
# Function to validate model-options provided by the user
# $1 must be a part of $2 list, else the script is stopped
# example: validateUserConfigValues "$queue" "$validateQueue"
function validateUserConfigValues {
validate="$(grep -w "$1" <<< "$2")"
if [ -z "${validate// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  echo "EXIT! > ""$1"" <  is not a valid value for parameter ""$3"", please change and restart prep_run_slm.sh"
  exit 1
 else
  if [ -z $1 ]
   then
    echo "EXIT!  The value for parameter ""$3"" is empty but must contain a valid string, please change and restart prep_run_slm.sh"
    exit 1
  fi
fi
return 0
}
# lists to validate against
validateForcingDataSet="fluxnet fluxnet3 cruncep face manip fert craine photacc spp1685 confin"
validateForcingMode="static transient"
validateCrunSites="all cue spp"
validateSlmModel="land plant soil canopy"
validateBNFScheme="none fixed standard optimal resistance unlimited"
validateMycScheme="no am ecm"
validateNDepositionScheme="none constant dynamic"
validatePDepositionScheme="none constant dynamic"
validateCanopyLayerScheme="standard fapar"
validateSBNLossScheme="none fixed dynamic"
validateBiomassAllocScheme="fixed dynamic optimal"
validateLeafStoichomScheme="fixed dynamic optimal"
validateVegDynamicsScheme="none population cohort"
validateSBModelScheme="simple_0d simple_1d jsm"
validateQueue="all io BTM mpi_large io_perform short"
validateEmailNotification="yes no"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validateOutputIntervalSpinup="none daily weekly monthly yearly"
validateSBAdsorpScheme="eca_full eca_part eca_none"


## set and create basic folders and files ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
mkdir -p "$bgiScratchUserFold"
slmBinary="../../x86_64-intel/bin/land.x"                  ## just one binary for all SLM models land plant soil canopy test_canopy
slmRunScriptJobFolder="./"
slmRunScriptJobFilename="run_slm.sh"
slmCleanupScriptFolder="./"
slmCleanupScriptFilename="cleanup_slm_hpc_run.sh"
slmNamelistFolder="../namelist/"
slmNamelistFilename="namelist.slm"
slmLctlibFolder="../../data/"
slmLctlibFilename="lctlib_slm_nlct8.def"
slmAdditionalDir="../additional_scripts/"
slmParamFigureFilename="create_parameter_summary_figs.R"    ## this filename is also used/defined in cleanup_slm_hpc_run.sh script
slmDefOutOverviewFilename="postprocessing_site_overview.R"  ## this filename is also used/defined in cleanup_slm_hpc_run.sh script
slmParamSensiFolder="../parameter_list/"
slmParamListRunFilename="parameter_slm_run.list"
slmParamTableFilename="parameter_table.txt"               ## this filename is also used/defined in cleanup_slm_hpc_run.sh script
# runscript cleanupscript namelist lctlib & parameter-variation files
 slmRunScriptJob="$slmRunScriptJobFolder""$slmRunScriptJobFilename"         ## runscript runing a job
 slmCleanupScript="$slmCleanupScriptFolder""$slmCleanupScriptFilename"      ## cleanup script started as extra job after simulations are done
 slmNamelist="$slmNamelistFolder""$slmNamelistFilename"                     ## just one namelist for all models
 slmLctlib="$slmLctlibFolder""$slmLctlibFilename"                           ## Landcovertype library for SLM with 8 PFTs
 slmParamListRun="$slmParamSensiFolder""$slmParamListRunFilename"           ## kind of namelist with only the parameter of which the value is varied
 slmParamTable="$slmParamSensiFolder""$slmParamTableFilename"               ## table with values of the parameter(s)
 slmParamFigureScript="$slmAdditionalDir""$slmParamFigureFilename"          ## R script to create overview figures from parameter variation simulations
 slmDefOutOverviewScript="$slmAdditionalDir""$slmDefOutOverviewFilename"    ## R script to create a .pdf with figs from one output file for all sites
# R scripts
 slmRscript="../r_script/postprocessing_""$slmModel"".R"                          ## R-Script for postprocessing output data
 slmRscript4fluxfile="../r_script/convert_slm2FLUXNET.R"                          ## R-Script for converting output for FLUXNET 
 slmRscript4fluxdiag="../r_script/generate_FLUXNET_diagnostics.R"                 ## R-Script for FLUXNET diagnostics 
 slmRscript4fluxBarPlots="../r_script/sim_fluxnet_diagnostics-process_timestep.R" ## R-Script for postprocessing FLUXNET diagnostics as barplots
 slmRscriptCompCruncepGfdb1="../r_script/sim_cruncep_vs_gfdb-process_daily.R"     ## R-Script for postprocessing simulation data for CRUNCEP vs GFDB
 slmRscriptCompCruncepGfdb2="../r_script/sim_cruncep_vs_gfdb-create_figs.R"       ## R-Script for creating comparision figures CRUNCEP vs GFDB
 slmRscriptCompCruncepGfdb3="../r_script/sim_cruncep_all_vs_gfdb-create_figs.R"   ## R-Script for creating comparision figures CRUNCEP vs GFDB

## experimental ID & simulation specific folders at scratch ##
# if experimental ID is empty, use date as: YYYYMMDD_HHMMSS
if [ -z $experimentalID ]
 then
  experimentalID=$(date +%Y%m%d_%H%M%S)
fi
# define job name
jobName="exp_""$experimentalID"                       # name of the job(array) -> must start with a character make work (see below): bsub -w "ended($jobName)"
simFold="$bgiScratchUserFold""$jobName""/"
simDefaultFiles="default_files/"
simParamSensiFiles="parameter_sensitivity_files/"     # this folder name is also used/defined in the cleanup_slm_hpc_run.sh script
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
      exit 1
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
      exit 1
  esac
fi


## validate user-defined model-options ##
validateUserConfigValues "$forcingDataSet" "$validateForcingDataSet" "forcingDataSet"
validateUserConfigValues "$forcingMode" "$validateForcingMode" "forcingMode"
validateUserConfigValues "$crunSites" "$validateCrunSites" "crunSites"
validateUserConfigValues "$customSiteList" "$validateFlags" "customSiteList"
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$queue" "$validateQueue" "queue"
validateUserConfigValues "$emailNotification" "$validateEmailNotification" "emailNotification"
validateUserConfigValues "$specialFluxnetOutput" "$validateFlags" "specialFluxnetOutput"
validateUserConfigValues "$bnf_scheme" "$validateBNFScheme" "bnf_scheme"
validateUserConfigValues "$flag_mycorrhiza" "$validateFlags" "flag_mycorrhiza"
validateUserConfigValues "$flag_mycorrhiza_org" "$validateFlags" "flag_mycorrhiza_org"
validateUserConfigValues "$flag_mycorrhiza_prim" "$validateFlags" "flag_mycorrhiza_prim"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$sb_model_scheme" "$validateSBModelScheme" "sd_model_scheme"
validateUserConfigValues "$sb_nloss_scheme" "$validateSBNLossScheme" "sb_nloss_scheme"
validateUserConfigValues "$veg_dynamics_scheme" "$validateVegDynamicsScheme" "veg_dynamics_scheme"
validateUserConfigValues "$flag_stand_harvest" "$validateFlags" "flag_stand_harvest"
validateUserConfigValues "$flag_sb_prescribe_nh4" "$validateFlags" "flag_sb_prescribe_nh4"
validateUserConfigValues "$flag_sb_prescribe_no3" "$validateFlags" "flag_sb_prescribe_no3"
validateUserConfigValues "$flag_sb_prescribe_po4" "$validateFlags" "flag_sb_prescribe_po4"
validateUserConfigValues "$flag_t_resp_acclimation" "$validateFlags" "flag_t_resp_acclimation"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_sb_jsm_transport" "$validateFlags" "flag_sb_jsm_transport"
validateUserConfigValues "$flag_sb_jsm_litter_input" "$validateFlags" "flag_sb_jsm_litter_input"
validateUserConfigValues "$flag_sb_jsm_OM_sorption" "$validateFlags" "flag_sb_jsm_OM_sorption"
validateUserConfigValues "$sb_adsorp_scheme" "$validateSBAdsorpScheme" "sb_adsorp_scheme"
validateUserConfigValues "$biomass_alloc_scheme" "$validateBiomassAllocScheme" "biomass_alloc_scheme"
validateUserConfigValues "$leaf_stoichom_scheme" "$validateLeafStoichomScheme" "leaf_stoichom_scheme"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$output_interval_pool_spinup" "$validateOutputIntervalSpinup" "output_interval_pool_spinup"
validateUserConfigValues "$output_interval_flux_spinup" "$validateOutputIntervalSpinup" "output_interval_flux_spinup"
validateUserConfigValues "$flag_transientSim_CO2const" "$validateFlags" "flag_transientSim_CO2const"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$flag_read_dC13" "$validateFlags" "flag_read_dC13"
validateUserConfigValues "$flag_read_DC14" "$validateFlags" "flag_read_DC14"
validateUserConfigValues "$useSiteSpecificAWC" "$validateFlags" "useSiteSpecificAWC"
validateUserConfigValues "$paramSensitivitySim" "$validateFlags" "paramSensitivitySim"
validateUserConfigValues "$paramChangesProportional" "$validateFlags" "paramChangesProportional"
validateUserConfigValues "$compareSimCRUNCEPvsGFDB" "$validateFlags" "compareSimCRUNCEPvsGFDB"



## set and create additional basic folders (within the basic simulation folder) ##
mkdir -p "$simFold"                           # "$bgiScratchUserFold""$jobName""/"  e.g. /Net/Groups/BGI/scratch/""$username""/slm/$jobName/
mkdir -p "$simFoldDefaultFiles"               # "$simFold""$simDefaultFiles"        e.g. /Net/Groups/BGI/scratch/""$username""/slm/$jobName/default_files/
mkdir -p "$simFoldOut"                        # "$simFold""out""/"                  e.g. /Net/Groups/BGI/scratch/""$username""/slm/$jobName/out/
mkdir -p "$simFoldErr"                        # "$simFold""error""/"                e.g. /Net/Groups/BGI/scratch/""$username""/slm/$jobName/error/
if [ $paramSensitivitySim == "T" ]; then
  mkdir -p "$simFoldParamSensiFiles"          # "$simFold""$simParamSensiFiles"     e.g. /Net/Groups/BGI/scratch/""$username""/slm/$jobName/parameter_sensitivity_files/
fi


## site data-set: FLUXNET / CRUNCEPv7 / FACE / MANIP / ... ##
case "$forcingDataSet" in
"fluxnet")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/"
  #siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET2/"
  siteDataObsFold="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/FLUXNET/"   # evaluation data; 'FLUXNET' is a symlink
  ;;
"fluxnet3")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET3/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"cruncep")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CRUNCEPv7/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"face")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FACE/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"manip")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/MANIP/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"fert")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/FERT/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"craine")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CRAINE/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"photacc")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/PHOTACC/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"spp1685")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/SPP_1685/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
"confin")
  siteDataFold="/Net/Groups/BGI/data/Model/SLM1/input/point/CONFIN/"
  siteDataObsFold="no such evaluation data"  # but this variable must contain something
  ;;
esac


## QUINCY model forcing: land/plant/canopy/soil (soil-forcing with 15 layers available)
# the soil-forcing dir is set below for 'case "$slmModel" in "soil")'
# met-forcing
metForcingDir="met_forcing/"
case "$slmModel" in
"land")
  siteDataFoldForcing="$siteDataFold""$metForcingDir"
  ;;
"plant")
  siteDataFoldForcing="$siteDataFold""$metForcingDir"
  ;;
"canopy")
  siteDataFoldForcing="$siteDataFold""$metForcingDir"
  ;;
"soil")
  # test number of soil layers against the default value 15; soil-forcing is available for 15 sl only
  if [ "$nsoil_layer" != 15 ]; then
    echo " "
    echo " You aim at using ""$nsoil_layer"" soil layers, where soil forcing is available for 15 sl only."
    echo " nsoil_layer is now changed to 15."
    nsoil_layer=15
  fi
  # test soil depth against the value necessary to run 15 sl, i.e., 9.5 meters
  if [ "$soil_depth" != 9.5 ]; then
    echo " You aim at using 15 soil layers, which needs a soil depth of 9.5 meters. But your soil depth is set to ""$soil_depth""."
    echo " soil_depth is now changed to 9.5 meters."
    soil_depth=9.5
  fi
  #
  # set soil-forcing directory
  #
  # static/transient + revNumber
  soilForcStatDirBase="soil_forcing_stat_sl15"
  soilForcTransDirBase="soil_forcing_trans_sl15"   # available only for a limited number of site-sets
  # cue/spp static forcing for CRUNCEPv7
  if [ "$forcingDataSet" == "cruncep" ]; then
    if [ "$crunSites" == "cue" ]; then
      soilForcStatDirBase="soil_forcing_stat_cue_sl15"
    elif [ "$crunSites" == "spp" ]; then
      soilForcStatDirBase="soil_forcing_stat_spp_sl15"
    fi
  fi
  # give static/transient value to a unified variable
  if [ "$forcingMode" == "static" ]; then
    soilForcDirBase="$soilForcStatDirBase"
    soilForcRevNum="$soilForcStatRevNum"
  else # transient
    soilForcDirBase="$soilForcTransDirBase"
    soilForcRevNum="$soilForcTransRevNum"
  fi
  # set soil forcing directory by adding the revision number
  #   current fallback for static soil-forcing: if no revision number is provided, a symlink points to the default soil-forcing
  if [ -z "$soilForcRevNum" ]; then   # if "$soilForcStatRevNum" == ""  (empty)
    if [ "$forcingMode" == "static" ]; then
      echo " The current default soil forcing is used. See symlink: BGI/data/Model/SLM1/input/point/[site-set]/""$soilForcStatDirBase"
      soilForcingDir="$soilForcDirBase""/"
    else # transient
      echo " EXIT! when using transient soil-forcing a revision number must be provided"
      exit 1
    fi
  else
    soilForcingDir="$soilForcDirBase""_r""$soilForcRevNum""/"
  fi
  siteDataFoldForcing="$siteDataFold""$soilForcingDir"
  # test if soil forcing is available for this site set
  soilForcingAvailable="T"
  if [ ! -d "$siteDataFoldForcing" ]; then
    soilForcingAvailable="F"    # directory not present
  else
    if [ -z "$(ls -A $siteDataFoldForcing)" ]; then
      soilForcingAvailable="F"  # directory empty
    fi
  fi
  # soil forcing availability feedback/action
  if [ "$soilForcingAvailable" == "T" ]; then
    # soil forcing available
    echo " Running the QUINCY soil model with forcing from:"
    echo " ""$siteDataFoldForcing"
  else
    echo " No soil-forcing available in the directory you provided:"
    echo " ""$siteDataFoldForcing"
    echo "EXIT!"
    exit 1
  fi
  ;;
esac


## test site list: default or custom ##
# parameter sensitivity: set customSiteList to "T" and delete any name in singleSiteName (see the IF[] for starting the job very below)
if [ $paramSensitivitySim == "T" ]; then
  echo ""
  echo " You aim at parameter sensitivity simulations for the site: ""$paramSensitivitySiteName"
  echo " Simulations will be done for each of the parameter combination you specified in contrib/parameter_list/parameter_table.txt"
  if [ $paramChangesProportional == "T" ]; then
    echo " Parameter values will be changed by a certain proportion [0<x<y], where y is not limited (i.e. you set paramChangesProportional=T)"
  else
    echo " Parameter values will be changed using absolute values as defined in parameter_table.txt (i.e. you set paramChangesProportional=F)"
  fi
  singleSiteName=""
fi
# test default and custom site list
if [ $customSiteList == "T" ]; then
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
    exit 1
   else 
    if [ ! -s "$siteListFileName" ] # if file is empty
     then 
    echo "EXIT! ""$siteListFileName"" has no content, i.e. is empty! You set the option customSiteList to T."
    exit 1
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
      # test if a forcing file for this site is available in the met_forcing / soil_forcing_sl15 directory
      # currently static & transient forcing is not tested individually  TODO
      if ! ls -la $siteDataFoldForcing | grep $siteName > /dev/null  # if forcing file not found
        then
          echo "EXIT! The site ""$siteName"" from your custom site_list.dat does not exist in the forcingDataSet folder of: ""$forcingDataSet"" !"
          rm -r "$simFold"  # deleting the already created simulation folder
          exit 1
      fi
  done
  # copy custom site-list and set link to file accordingly
  cp "$siteListFileName" "$simFoldDefaultFiles"
  siteList="$simFoldDefaultFiles""$siteListFileName"
# the default case: $customSiteList == "F"
else 
  ## select the correct site information file (a.k.a. "all_site_list.dat")
  ## file names differ across forcingDataSet & forcingMode
  case "$forcingDataSet" in
  "fluxnet")
    if [ $forcingMode == "static" ]; then
      siteListFileName="flux_all_sites_1995-2006_list.dat"
    else 
      siteListFileName="flux_all_sites_1901-2006_list.dat"
    fi
    ## SLM baseline paper 2018 using a specific set of fluxnet sites per default for usecases 80 81 82
    if [ "$usecase" != "none" ]; then  # if a $usecase is defined
      if [ "$usecase" -gt 79 ] && [ "$usecase" -lt 83 ]; then  # if usecase 80 or 81 or 82
        siteListFileName="flux_baseline_sites_1901-2006_list.dat"
        echo " You use u80/81/82. Setting site list file to flux_baseline_sites_1901-2006_list.dat."
      fi
    fi
    ;;
  "fluxnet3")
    if [ $forcingMode == "static" ]; then
      siteListFileName="flux3_all_sites_1991-2014_list.dat"
    else 
      siteListFileName="flux3_all_sites_1901-2014_list.dat"
    fi
    ;;
  "cruncep")
    case "$crunSites" in
    "all")
      if [ $forcingMode == "static" ]; then
        siteListFileName="crun_all_gfdb_sites_1984-2013_list.dat"
      else
        siteListFileName="crun_all_gfdb_sites_1901-2015_list.dat"
      fi
      # no soil forcing for "cruncep all" yet
      ;;
    "cue")
      if [ $forcingMode == "static" ]; then
        siteListFileName="crun_cue_gfdb_sites_1984-2013_list.dat"
      else
        siteListFileName="crun_cue_gfdb_sites_1901-2015_list.dat"
      fi
      ;;
    "spp")
      if [ $forcingMode == "static" ]
       then
        siteListFileName="crun_spp_sites_1984-2013_list.dat"
      else
        siteListFileName="crun_spp_sites_1901-2015_list.dat"
      fi
      ;;
    esac
    ;;
  "face")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="face_sites_1996-2008_list.dat"
    else
      siteListFileName="face_sites_1901-2015_list.dat"
    fi
    ;;
  "manip")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="manip_sites_2009-2019_list.dat"
    else
      echo " Sorry. Transient forcing for MANIP sites not yet available."
      echo "EXIT!"
      exit 1
      #siteListFileName="manip_sites_xxxx-2017_list.dat"
    fi
    ;;
  "fert")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="fert_sites_1901-2015_list.dat"
      echo " "
      echo " Note: you are running a static simulation with (a) FERT site(s) using (transient) forcing 1901-2015 by default."
      echo " "
    else
      siteListFileName="fert_sites_1901-2015_list.dat"
    fi
    ;;
  "craine")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="craine_sites_1901-2015.dat"
      echo " "
      echo " Note: you are running a static simulation with (a) CRAINE site(s) using (transient) forcing 1901-2015 by default."
      echo " "
    else
      siteListFileName="craine_sites_1901-2015.dat"
    fi
    ;;
  "photacc")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="photacc_sites_1901-2015.dat"
      echo " "
      echo " Note: you are running a static simulation with (a) PHOTACC site(s) using (transient) forcing 1901-2015 by default."
      echo " "
    else
      siteListFileName="photacc_sites_1901-2015.dat"
    fi
    ;;
  "spp1685")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="spp_1685_sites_1984-2013_list.dat"
    else 
      siteListFileName="spp_1685_sites_1901-2015_list.dat"
    fi
    ;;
  "confin")
    if [ $forcingMode == "static" ]
     then
      siteListFileName="confin_sites_1997-2014_list.dat"
    else 
      siteListFileName="confin_sites_1901-2014_list.dat"
    fi
    ;;
  esac
  siteList="$siteDataFold""$siteListFileName"
fi


## apply setting for comparing simulated model results of CRUNCEP sites with Global Forest Data Base ##
# this may overwrite some options specifies above !
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then
  if [ "$forcingDataSet" != "cruncep" ]; then
    echo "EXIT!  You want to run a compareSimCRUNCEPvsGFDB simulation with model-data comparison, but the site set is not cruncep."
    exit 1
  fi
  echo " You are running the CRUNCEP vs GFDB simulation for the cruncep > ""$crunSites"" < site selection."
  # set path/name for files with field data, for checking their presence below, as long as dailyOutput="T"
  gfdbFileCUE="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_cue_sites_01082017.csv"
  gfdbFileALL="/Net/Groups/BGI/data/Model/SLM1/evaluation/point/GFDB/global_forest_database_all_sites_01082017.csv"
  dailyOutput="T"
  if [ "$output_interval_pool" != "daily" ]; then 
    echo " You are running the CRUNCEP vs GFDB simulation with ""$output_interval_pool"" pool output."
    echo " The sim_cruncep_vs_gfdb-process_daily.R needs daily data - hence, the automatic model-data comparison may not work."
    dailyOutput="F"
  fi
  if [ "$output_interval_flux" != "daily" ]; then 
    echo " You are running the CRUNCEP vs GFDB simulation with ""$output_interval_flux"" flux output."
    echo " The sim_cruncep_vs_gfdb-process_daily.R needs daily data - hence, the automatic model-data comparison may not work."
    dailyOutput="F"
  fi
  # disable parameter sensitivity runs
  paramSensitivitySim="F"
  # disable special fluxnet-output
  specialFluxnetOutput="F"
  if [ $forcingMode == "static" ]
   then
    # output starts right after the spinup has finished (output is used for simulation vs GFDB comparision)
    let simulation_length_number=spinupLengthGFDB+simulation_length_number
    let output_start_first_day_year=spinupLengthGFDB+output_start_first_day_year
    let output_end_last_day_year=spinupLengthGFDB+output_end_last_day_year
    if [ "$veg_dynamics_scheme" == "cohort" ] && [ "$spinupLengthGFDB" -gt 0 ]; then
      echo " Note: the veg_dynamics_scheme cohort does not work with the static forcing mode + spinup. Sorry."
      echo " Either set spinupLengthGFDB to zero or use the transient mode for spinup."
      echo "EXIT!"
      exit 1
    fi
  else
    # transient simulations:
    # overwrites the above setting for transient_spinup_years
    transient_spinup_years="$spinupLengthGFDB"
    # simulation_length_number does not change, its the number of years after the transient_spinup_years
    # neither do the output_start* & output_end* years change
    # note, in transient simulations the transient_spinup_years are counted within SLM from negative value to zero, which is when the simulation_length_number starts
  fi
  # checking GFDB field-data files availability if dailyOutput="T" (if no daily output is produced, the automatic model-data comparison may not work anyway)
  if [ "$dailyOutput" == "T" ]; then
    if [ ! -e "$gfdbFileCUE" ] || [ ! -e "$gfdbFileALL" ]; then  # if one of the files does not exists
      echo "EXIT! The GFDB field-data files are not available! These are needed by the sim_cruncep*_vs_gfdb-create_figs.R scripts."
      echo " File locations should be: ""$gfdbFileCUE"" and ""$gfdbFileALL"
      exit 1
    else
      echo " "
      echo " You started simulations to compare CRUNCEP with GFDB data using daily output and a spin-up time of ""$spinupLengthGFDB"" years."
      if [ $forcingMode == "static" ]; then
        echo " The comparison will be based on daily SLM output of "$(($output_end_last_day_year - $output_start_first_day_year + 1))" years after the spin-up (static forcing mode)."
      fi
    fi
  fi # check "$dailyOutput" == "T"
fi


## evaluate options: forcingMode == transient ##
#  some settings of the transient simulations are affected by the compareSimCRUNCEPvsGFDB simulation settings if compareSimCRUNCEPvsGFDB == "T"
if [ $forcingMode == "transient" ] && [ "$slmModel" != "soil" ]; then
  echo " "
  echo " You use QUINCY ""$slmModel"" and set the forcingMode to transient, with a spinup time of ""$transient_spinup_years"" years that uses forcing input from year ""$transient_spinup_start_year"" to ""$transient_spinup_end_year""."
  echo " The subsequent simulation may cover ""$simulation_length_number"" years, starting with forcing data from year ""$transient_simulation_start_year""."

  # spinup length
  if [ $transient_spinup_years -lt 1 ]; then
    echo " Ohh, you set the spinup time to a value less than 1 year, which is not recommended. Re-setting it to 1 now."
    transient_spinup_years=1
  fi
  if [ $transient_spinup_years -lt 50 ]; then
    echo " Warning: You set the spinup length to ""$transient_spinup_years"" years, which might be to low for sufficient results (default would be 500 years)."
  fi
  # starting year 
  if [ $transient_spinup_start_year != "1901" ]; then
    echo " You set the transient_spinup_start_year to ""$transient_spinup_start_year"", which is currently not supported. Sorry."
    echo " Setting transient_spinup_start_year=1901."
    transient_spinup_start_year=1901
  fi
  if [ $transient_simulation_start_year != "1901" ]; then
    echo " You set the transient_simulation_start_year to ""$transient_simulation_start_year"", which is currently not supported. Sorry."
    echo " Setting transient_simulation_start_year=1901."
    transient_simulation_start_year=1901
  fi
  # length of forcing period repeatedly used for spinup
  transient_test_spinup_forcing_period=$((transient_spinup_end_year - transient_spinup_start_year + 1))
  if [ $transient_test_spinup_forcing_period -lt 2 ] || [ $transient_test_spinup_forcing_period -gt 49 ]; then
    echo " You set the length of forcing period repeatedly used for spinup to ""$transient_test_spinup_forcing_period""."
    echo " It is recommended to use a period length (transient_spinup_end_year - transient_spinup_start_year) between 2 and 50 years."
    transient_spinup_end_year=$((transient_spinup_start_year + 29))
    echo " The new (corrected) value of transient_spinup_end_year is ""$transient_spinup_end_year"" setting the period length to the default of 30 years."
  fi
# transient & QUINCY soil model
elif [ $forcingMode == "transient" ] && [ "$slmModel" == "soil" ]; then
  echo " "
  echo " You use QUINCY ""$slmModel"" and set the forcingMode to transient, with a spinup time of ""$transient_spinup_years"" years."
  echo " The transient simulation (after) spinup) may cover one full forcing period, which is a site-specific parameter (fluxnet/3 site-sets )."  
  echo " EXIT! not yet implemented"
  exit 1
# TODO  add code for soil transient simulations here, and to the run_slm.sh !
# forcingMode == "static"
else  
  # set transient* parameters to default values, to avoid potentially quirky values here, but do not want do do checking either
  transient_spinup_start_year=1901
  transient_spinup_end_year=1930
  transient_spinup_years=500
  transient_simulation_start_year=1901
fi


## transient forcing and cruncep vs gfdb
if [ $compareSimCRUNCEPvsGFDB == "T" ] &&  [ $forcingMode == "transient" ]; then
  let startYearTransOut=transient_simulation_start_year+output_start_first_day_year-1
  echo " The comparison will be based on daily SLM output of "$(($output_end_last_day_year - $output_start_first_day_year + 1))" years starting in year ""$startYearTransOut"" (transient forcing mode)."  
fi  
  

## transient forcing and fluxnet simulation with 
# note: actual simulation length 1901 to x varies with the site, because the available site data does start at diff years and has different lengths
#       this is dealt with in run_slm.sh
fluxnetTransQLand="F"
if [ $forcingMode == "transient" ] && [ "$slmModel" != "soil" ]; then
  if [ $forcingDataSet == "fluxnet" ] || [ $forcingDataSet == "fluxnet3" ]; then
    echo " You started a transient simulation with QUINCY ""$slmModel"" and fluxnet forcing: ""$forcingDataSet"
    echo " Note: the simulation (starting in 1901 after the spinup) ends with the last year of data available for the particular site."
    echo " Hence, the actual simulation length may vary with the fluxnet sites."
    echo " You set the simulation length to: ""$simulation_length_number"" starting in 1901. This may change for most sites."
    fluxnetTransQLand="T"
  fi
fi


## an identifier (additional to siteName) for static and transient forcing files -> passed to run_slm.sh (via config file) ##
case "$forcingDataSet" in
"fluxnet")
  if [ $forcingMode == "static" ]; then
    forcingModeFileID="_s_"
  else  # valid for land & soil
    forcingModeFileID="_t_"
  fi
  # QUINCY soil model: no specific forcing-file indentifier needed; static & transient soil-forcing are stored in individual dirs
  if [ "$slmModel" == "soil" ]; then forcingModeFileID="_"; fi
  ;;
"fluxnet3")
  if [ $forcingMode == "static" ]; then
    forcingModeFileID="_s_"
  else  # valid for land & soil
    forcingModeFileID="_t_"
  fi
  # QUINCY soil model: no specific forcing-file indentifier needed; static & transient soil-forcing are stored in individual dirs
  if [ "$slmModel" == "soil" ]; then forcingModeFileID="_"; fi
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
"manip")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_2009"         # ES-MaG_2009-2017.dat & ES-MaT_2009-2017.dat
  else
    echo " ... set forcing filename ID for manip transient"
    echo " EXIT1"
    exit 1
    forcingModeFileID="_xxxx"
  fi
  ;;
"fert")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_1901"
  else
    forcingModeFileID="_1901"
  fi
  ;;
"craine")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_1901"
  else
    forcingModeFileID="_1901"
  fi
  ;;
"photacc")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_1901"
  else
    forcingModeFileID="_1901"
  fi
  ;;
"spp1685")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_1984"
  else
    forcingModeFileID="_1901"
  fi
  ;;
"confin")
  if [ $forcingMode == "static" ]
   then
    forcingModeFileID="_s_"
  else
    forcingModeFileID="_t_"
  fi
  ;;
esac


## validate experimental ID for alphanumeric only ##
if [ ! -z $(echo $jobName | tr -d "[:alnum:]" | tr -d "_" | tr -d "-") ]  # if removing alphanum _ - from jobName returns an non-empty string" then exit
 then
  echo "EXIT! Please use an experimental ID that is alphanumeric, i.e. without any special characters (_ and - are ok though)."
  exit 1
fi


## test number of soil layers against min and max values (5:20)
if [ "$nsoil_layer" -lt 5 ]
 then
  echo " "
  echo " You aim at using ""$nsoil_layer"" soil layers, where the minimum number is 5. nsoil_layer is now set to 5"
  nsoil_layer=5
 else
  if [ "$nsoil_layer" -gt 20 ]
   then
    echo " "
    echo " You aim at using ""$nsoil_layer"" soil layers, where the maximum number is 20. nsoil_layer is now set to 20"
    nsoil_layer=20
  fi
fi


## test cohort-harvesting interval & options related to veg_dynamics_scheme="cohort"
if [ "$veg_dynamics_scheme" == "cohort" ]
 then
   # cohort_harvest_interval range
   if [ "$cohort_harvest_interval" -lt 20 ] || [ "$cohort_harvest_interval" -gt 140 ]; then
     echo " Note. The cohort-harvesting interval might be set to a value from 20 to 140. You set it to: ""$cohort_harvest_interval""."
   fi
   # cohort_harvest_interval minimum
   if [ "$cohort_harvest_interval" -lt 1 ]; then
     echo " The cohort-harvesting interval is set to a value smaller than 1."
     echo "EXIT!"
     exit 1
   fi
   # stand harvest may be set to F
   if [ "$flag_stand_harvest" == "T" ]; then
     echo " The stand harvest may not be enabled when using the veg_dynamics_scheme cohort. Setting flag_stand_harvest to F now."
     echo " This is because the cohort is harvested automatically at the interval of cohort_harvest_interval."
     flag_stand_harvest="F"
   fi
fi


## set mycorrhiza flag_myco_* flags to FALSE if mycorrhiza are disabled overall
if [ $flag_mycorrhiza == "F" ]
 then
  flag_mycorrhiza_org="F"
  flag_mycorrhiza_prim="F"
fi


## validate options linked with specialFluxnetOutput a.k.a. fluxnet diagnostics simulations ##
#  requires transient forcing from 1901 to ~2015 !
#  in addition to the normal transient simulation:
#    the last full forcing period is written out at timestep 
#    extra R scripts are started per site after the simulation comparing simulation results with empirical data
if [ $specialFluxnetOutput == "T" ]; then
  # check site data set
  if [ $forcingDataSet == "fluxnet" ]; then
    echo " Running fluxnet diagnostics simulation with site-set: FLUXNET2"
  elif [ $forcingDataSet == "fluxnet3" ]; then
    echo " Running fluxnet diagnostics simulation with site-set: FLUXNET3"
  else
    echo "EXIT! You might want to use FLUXNET/3 sites when running fluxnet diagnostics simulation :) ."
    exit 1
  fi
  # check QUINCY model (fluxnet/3 site-set already tested above)
  if [ "$slmModel" == "soil" ]; then
    echo " You are trying to run the QUINCY soil model with a fluxnet diagnostics simulation setup (uses transient forcing)."
    echo "EXIT! fluxnet diagnostics simulationa require transient forcing from 1901 to ~2015. Such soil_forcing is not available for now."
    rm -r "$simFold"  # deleting the already created simulation folder
    exit 1
  fi
  if [ "$slmModel" != "land" ]; then
    echo "EXIT! The fluxnet diagnostics simulations should only run with the QUINCY land model, not the ""$slmModel"" model."
    rm -r "$simFold"  # deleting the already created simulation folder
    exit 1
  fi
  # message to screen
  echo " "
  echo " You just started simulations to generate FLUXNET/3 diagnostics from the last complete forcing period before the end of the ""$simulation_length_number"" years."
  # docu: if e.g. the simulation length is 80 years, the last full forcing period, considering 3 years of data, could be from year 76 to 78
  # if $singleSiteName is not empty
  if [ ! -z "$singleSiteName" ]; then
    echo " Simulation is done for one site: ""$singleSiteName"
    echo " "
  fi
  # parameter sensitivity simulations
  if [ $paramSensitivitySim == "T" ]; then
    echo " Multiple parameter sensitivity simulations are done for the site: ""$paramSensitivitySiteName"
    echo " "
  fi
  # $singleSiteName is empty and no parameter sensitivity simulations
  if [ -z "$singleSiteName" ] && [ $paramSensitivitySim == "F" ]; then
     echo " Simulations are done for all sites listed in: ""$siteListFileName"
     echo " "
  fi
fi


## validate allocated memory ##
lowestMem=200      # 200 MB
largestMem=30000   # 30  GB
# in case of using the "extrafast" option, the max memory allowed is increased, because that's how the "extrafast" option works
if [ "$efastOpt" == "TRUE" ]
 then
   largestMem=59500   # 59.5 GB | a node of the all queue has 60 GB memory installed
fi
# test memory empty or containing characters other then digits
if [ -z $memory ]  # returns TRUE if unset or empty
 then
  echo "EXIT! No memory defined. Please define the memory per job you think is useful (and within a reasonable range)."
  exit 1  
 else
  if [ ! -z $(echo $memory | tr -d "[:digit:]") ] # if "removing all digits" is not empty, i.e. other characters are included
   then
    echo "EXIT! The memory parameter contains characters other than digits. This is not working. Please use just numbers."
    exit 1
  fi
fi
# test against useful min and max values
if [ "$memory" -lt "$lowestMem" ]
 then
  echo "EXIT! Please increase memory, you allocated less than 200MB. Note, memory is defined with unit MegaByte [MB]."
  exit 1
 else
  if [ "$memory" -gt "$largestMem" ]
   then
    echo "EXIT! Please decrease memory, you allocated more than 30GB. Note, memory is allocated per job."
    exit 1
  fi
fi


## compile the bsub -R argument (Reserved) for "reserved memory" in MB
reservedArg="rusage[mem="$memory"]"
# in case the io_perform "virtual" queue is selected, add also processor information to the -R argument and change queue to 'io'
if [ "$queue" == "io_perform" ]
 then
  reservedArg="select[model==INTEL3] rusage[mem="$memory"]"
  queue="io"      # because there is actually no queue 'io_perform'
fi


## calculate memory (bsub -M argument) from MB -> KB as needed by OpenLava
memory=$((memory * 1000))


## validate allocated number of CPUs ##
lowestCpu=5      # 5   parallel jobs
largestCpu=350   # 350 parallel jobs
# test numberCpuMax empty or containing characters other then digits
if [ -z $numberCpuMax ]  # returns TRUE if unset or empty
 then
  echo "EXIT! No numberCpuMax values defined. Please define the number of CPUs that may run in parallel at the BGC Cluster."
  exit 1  
 else
  if [ ! -z $(echo $numberCpuMax | tr -d "[:digit:]") ] # if "removing all digits" is not empty, i.e. other characters are included
   then
    echo "EXIT! The numberCpuMax parameter contains characters other than digits. This is not working. Please use just numbers."
    exit 1
  fi
fi
# test against useful min and max values
if [ "$numberCpuMax" -lt "$lowestCpu" ]
 then
  echo "Note: This is a High Performance Cluster. You may consider running more than 4 simulations in parallel."
  echo " Starting your jobs anyway ... "
 else
  if [ "$numberCpuMax" -gt "$largestCpu" ]
   then
    echo "EXIT! Please decrease the number of parallel jobs, you aim at using more than 350 CPUs in parallel."
    echo "(Of course, technically you could use more than 350 CPUs. This limit is simply set to be nice with your colleagues.)"
    exit 1
  fi
fi

## convert $emailNotification into a bsub option ##
emailSwitch=""
if [ "$emailNotification" == "yes" ]
  then
    emailSwitch="-N"
fi

## test if namelist binary siteList R-Script (and parameter_slm_run.list + parameter_table.txt) exists ##
if [ ! -e "$slmBinary" ]     # exit if non-existing
  then
	echo "EXIT! binary is missing (the ifort compiled version is expected)"
	exit 1
fi
if [ ! -e "$slmNamelist" ]   # exit if non-existing
  then
	echo "EXIT! namelist is missing"
	exit 1
fi
if [ ! -e "$slmLctlib" ]   # exit if non-existing
  then
	echo "EXIT! lctlib is missing"
	exit 1
fi
if [ ! -e "$siteList" ]      # exit if non-existing
  then
	echo "EXIT! site-list-data file is missing"
	exit 1
fi	
if [ ! -e "$slmRscript" ]    # exit if non-existing
  then
	echo "EXIT! the R Script for postprocessing output data is missing"
	exit 1
fi
# in case of parameter sensitivity
if [ $paramSensitivitySim == "T" ]; then
  # parameter_slm_run.list - exit if non-existing
  if [ ! -e "$slmParamListRun" ]; then
  	echo "EXIT! the default parameter List (parameter_slm_run.list) is missing though needed for Parameter Sensitivity Simulations."
    echo " You may run the contrib/parameter_list/generate_parameter_table.sh to create this file."
  	exit 1
  fi
  # parameter_table.txt - exit if non-existing
  if [ ! -e "$slmParamTable" ]; then
  	echo "EXIT! the parameter Table (parameter_table.txt) is missing though needed for Parameter Sensitivity Simulations."
    echo " You may run the contrib/parameter_list/generate_parameter_table.sh to create the table."
  	exit 1
  fi	
fi

## validate $singleSiteName for running SLM at one specific site only ##
# validate if $singleSiteName is not empty
if [ ! -z "$singleSiteName" ]
 then
  allSites=$(cat "$siteList" | awk '{print $1}')
  validateUserConfigValues "$singleSiteName" "$allSites" "singleSiteName"
fi

## parameter sensitivity: validate $paramSensitivitySiteName and get the number of parameter-sets/simulations ##
if [ $paramSensitivitySim == "T" ]; then
  # validate $paramSensitivitySiteName
  allSites=$(cat "$siteList" | awk '{print $1}')
  validateUserConfigValues "$paramSensitivitySiteName" "$allSites" "paramSensitivitySiteName"
  # get line number of site in the siteList, needed for the run_slm.sh to get the site-specific data from the site list info file
  paramSensitivitySiteLineNumber=$(grep -n "$paramSensitivitySiteName" "$siteList" | cut -f1 -d ":")
  # calculate number of parameter-sets/simulations, needed also for the bsub job array
  paramSensi_NumberCombinations=$(wc -l "$slmParamTable" | awk '{print $1}')
  paramSensi_NumberCombinations=$(($paramSensi_NumberCombinations -1))  # not counting the header of the table
  # validate number of parameter combinations
  if [ $paramSensi_NumberCombinations -gt 1000 ]; then
    echo "EXIT !  You want to do ""$paramSensi_NumberCombinations"" simulations."
    echo "        This exceeds the max jobarray size. A number lower than 1001 is necessary."
    exit 1
   else
    echo " The number of combinations of parameter values (i.e. simulations) is: ""$paramSensi_NumberCombinations"
  fi
fi

############################################################################################################################
## copy namelist binary R-Scripts to simulation folder: /Net/Groups/BGI/scratch/""$username""/slm/$jobName/default_files/ ##
############################################################################################################################
cp "$slmRunScriptJob"           "$simFoldDefaultFiles"
cp "$slmCleanupScript"          "$simFoldDefaultFiles"
cp "$slmBinary"                 "$simFoldDefaultFiles"
cp "$slmNamelist"               "$simFoldDefaultFiles"
cp "$slmLctlib"                 "$simFoldDefaultFiles"
cp "$slmRscript"                "$simFoldDefaultFiles"
cp "$slmDefOutOverviewScript"   "$simFoldDefaultFiles"
if [ $specialFluxnetOutput == "T" ]; then 
  cp "$slmRscript4fluxfile"       "$simFoldDefaultFiles"
  cp "$slmRscript4fluxdiag"       "$simFoldDefaultFiles"
  cp "$slmRscript4fluxBarPlots"   "$simFoldDefaultFiles"
fi
if [ $paramSensitivitySim == "T" ]; then
  cp "$slmParamListRun"       "$simFoldParamSensiFiles"
  cp "$slmParamTable"         "$simFoldParamSensiFiles"
  cp "$slmParamFigureScript"  "$simFoldParamSensiFiles"
fi
if [ $compareSimCRUNCEPvsGFDB == "T" ]; then
  cp "$slmRscriptCompCruncepGfdb1" "$simFoldDefaultFiles"
  cp "$slmRscriptCompCruncepGfdb2" "$simFoldDefaultFiles"
  cp "$slmRscriptCompCruncepGfdb3" "$simFoldDefaultFiles"
fi

## syncronize all data with the file system ##
sleep 1 && sync

#############################################################################################
## define a function that manipulates the namelist && change namelist according to options ##
#############################################################################################
# model: slm_model_name = land plant soil canopy
# 2 do: add other options like optimization

 # this function is copied from run_slm.sh without any modification
 # function: pass the variable_name ($1) & the new_value ($2)
 # $slmNamelist refers to the namelist file name, assuming the namelist file is located in the current folder !
function changeNamelistValue {
oldValue=$(grep -w $1 "$simFoldDefaultFiles""$slmNamelistFilename" | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 "$simFoldDefaultFiles""$slmNamelistFilename" | cut -d : -f 1)                       # get line number of variable $1
sed -i "$lineNum"s/$oldValue/$2/ "$simFoldDefaultFiles""$slmNamelistFilename"                             # replace the old value in line $lineNum with $2
# explicitly testing and error catching to avoid using sync
sleep 0.1
newValue=$(grep -w $1 "$simFoldDefaultFiles""$slmNamelistFilename" | awk -F "=" '{print $2}' | tr -d '"') # get the new value of this variable
if [ "$2" != "$newValue" ]; then
  sleep 2
  newValue=$(grep -w $1 "$simFoldDefaultFiles""$slmNamelistFilename" | awk -F "=" '{print $2}' | tr -d '"') # get the new value of this variable
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
 #changeNamelistValue soil_theta_prescribe 1.


##################################################################################
## namelist manipulation                                                        ##
##################################################################################
# change values in the namelist copied to the simulation folder at the scratch drive
# this namelist is used for all sites
changeNamelistValue slm_model_name "$slmModel"
changeNamelistValue forcing_mode "$forcingMode"
changeNamelistValue bnf_scheme "$bnf_scheme" 
changeNamelistValue n_deposition_scheme "$n_deposition_scheme"
changeNamelistValue p_deposition_scheme "$p_deposition_scheme"
changeNamelistValue ncanopy "$ncanopy"
changeNamelistValue canopy_layer_scheme "$canopy_layer_scheme" 
changeNamelistValue sb_model_scheme "$sb_model_scheme" 
changeNamelistValue sb_nloss_scheme "$sb_nloss_scheme"
changeNamelistValue sb_bnf_scheme "$bnf_scheme" 
changeNamelistValue veg_dynamics_scheme "$veg_dynamics_scheme"
changeNamelistValue biomass_alloc_scheme "$biomass_alloc_scheme"
changeNamelistValue leaf_stoichom_scheme "$leaf_stoichom_scheme"
changeNamelistValue cohort_harvest_interval "$cohort_harvest_interval"
changeNamelistValue sb_adsorp_scheme "$sb_adsorp_scheme"
if [ "$forcingDataSet" == "spp1685" ]
  then
    changeNamelistValue flag_spp1685 ".TRUE."
  else
    changeNamelistValue flag_spp1685 ".FALSE."
fi
if [ $flag_mycorrhiza == T ]
  then
    changeNamelistValue flag_mycorrhiza ".TRUE."
  else
    changeNamelistValue flag_mycorrhiza ".FALSE."
fi
if [ $flag_mycorrhiza_org == T ]
  then
    changeNamelistValue flag_mycorrhiza_org ".TRUE."
  else
    changeNamelistValue flag_mycorrhiza_org ".FALSE."
fi
if [ $flag_mycorrhiza_prim == T ]
  then
    changeNamelistValue flag_mycorrhiza_prim ".TRUE."
  else
    changeNamelistValue flag_mycorrhiza_prim ".FALSE."
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
if [ "$paramChangesProportional" == "T" ]
  then
    changeNamelistValue set_parameter_values_proportional ".TRUE."
  else
    changeNamelistValue set_parameter_values_proportional ".FALSE."
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
if [ $flag_read_DC14 == T ]
  then
    changeNamelistValue flag_read_DC14 ".TRUE."
  else
    changeNamelistValue flag_read_DC14 ".FALSE."
fi
if [ $flag_sb_jsm_transport == T ]
  then
    changeNamelistValue flag_sb_jsm_transport ".TRUE."
  else
    changeNamelistValue flag_sb_jsm_transport ".FALSE."
fi
if [ $flag_sb_jsm_litter_input == T ]
  then
    changeNamelistValue flag_sb_jsm_litter_input ".TRUE."
  else
    changeNamelistValue flag_sb_jsm_litter_input ".FALSE."
fi
if [ $flag_sb_jsm_OM_sorption == T ]
  then
    changeNamelistValue flag_sb_jsm_OM_sorption ".TRUE."
  else
    changeNamelistValue flag_sb_jsm_OM_sorption ".FALSE."
fi
if [ $flag_transientSim_CO2const == T ]
  then
    changeNamelistValue flag_transient_sim_co2_const ".TRUE."
  else
    changeNamelistValue flag_transient_sim_co2_const ".FALSE."
fi
if [ $flag_oneYearForcing == T ]
  then
    changeNamelistValue flag_climate_forcing_one_year ".TRUE."
  else
    changeNamelistValue flag_climate_forcing_one_year ".FALSE."
fi
changeNamelistValue simulation_length_number "$simulation_length_number"
changeNamelistValue output_interval_pool "$output_interval_pool"
changeNamelistValue output_interval_flux "$output_interval_flux"
changeNamelistValue output_interval_pool_spinup "$output_interval_pool_spinup"
changeNamelistValue output_interval_flux_spinup "$output_interval_flux_spinup"
changeNamelistValue output_start_first_day_year "$output_start_first_day_year"
changeNamelistValue output_end_last_day_year "$output_end_last_day_year"
changeNamelistValue soil_awc_prescribe "$soil_awc_prescribe"
changeNamelistValue soil_theta_prescribe "$soil_theta_prescribe"
changeNamelistValue soil_depth "$soil_depth"
changeNamelistValue nsoil_energy "$nsoil_layer"
changeNamelistValue nsoil_water "$nsoil_layer"
changeNamelistValue transient_spinup_start_year "$transient_spinup_start_year"
changeNamelistValue transient_spinup_end_year "$transient_spinup_end_year"
changeNamelistValue transient_spinup_years "$transient_spinup_years"
changeNamelistValue transient_simulation_start_year "$transient_simulation_start_year"
## namelist options for soil forcing creation (BASE_CTL):
## u01: static, one year, all site-sets | u08: transient, fluxnet-forcing years, fluxnet/3 sites
createSoilForcFlag="F"    # default: F | needed for u08 in run_slm.sh
if [ $usecase == "01" ]; then                                     # usecase file: 01_generate_soil_forcing_def
  changeNamelistValue flag_generate_soil_forcing_stat ".TRUE."    # default is .FALSE.
  createSoilForcFlag="T"                                          # needed only for u08, but set T here anyway
fi
if [ $usecase == "08" ]; then                                     # usecase file: 08_generate_flux_trans_soil_forcing_def
  changeNamelistValue flag_generate_soil_forcing_trans ".TRUE."   # default is .FALSE.
  createSoilForcFlag="T"                                          # disable all "fluxnet diagnostics" scripts & stuff in run_slm.sh
fi  

#########################################################################################################
## pass config info and filenames to run_slm.sh via an extra file ("$simFoldDefaultFiles""config.txt") ##
#########################################################################################################
# /Net/Groups/BGI/scratch/""$username""/slm/$jobName/default_files/config.txt
echo "$siteList"                                                  >  "$simFoldDefaultFiles""config.txt"   # 1
echo "$siteDataFoldForcing"                                       >> "$simFoldDefaultFiles""config.txt"
echo "$slmModel"                                                  >> "$simFoldDefaultFiles""config.txt"
echo "$forcingDataSet"                                            >> "$simFoldDefaultFiles""config.txt"
echo "$forcingModeFileID"                                         >> "$simFoldDefaultFiles""config.txt"   # 5
echo "$specialFluxnetOutput"                                      >> "$simFoldDefaultFiles""config.txt"
echo "$siteDataObsFold"                                           >> "$simFoldDefaultFiles""config.txt"
echo "$useSiteSpecificAWC"                                        >> "$simFoldDefaultFiles""config.txt"
echo "$paramSensitivitySim"" ""$simParamSensiFiles"               >> "$simFoldDefaultFiles""config.txt"
echo "$compareSimCRUNCEPvsGFDB"                                   >> "$simFoldDefaultFiles""config.txt"   # 10
echo "$fluxnetTransQLand"" ""$transient_simulation_start_year"    >> "$simFoldDefaultFiles""config.txt"
echo "$createSoilForcFlag"                                        >> "$simFoldDefaultFiles""config.txt"
#echo "$"                                                         >> "$simFoldDefaultFiles""config.txt"
# pass extra info about parameter sensitivity simulations to folder
if [ $paramSensitivitySim == "T" ]; then
  echo "$paramSensitivitySiteLineNumber"         >  "$simFoldParamSensiFiles""param_config.txt"
  echo "$slmParamListRunFilename"               >>  "$simFoldParamSensiFiles""param_config.txt"
  echo "$slmParamTableFilename"                 >>  "$simFoldParamSensiFiles""param_config.txt"
  #echo "$"                                     >>  "$simFoldParamSensiFiles""param_config.txt"
  #echo "$"                                     >>  "$simFoldParamSensiFiles""param_config.txt"
fi

## syncronize all data with the file system ##
sleep 1 && sync

##################################################################################
## start job with simulations by starting the job(array) -> run_slm.sh          ##
##################################################################################
# run_slm.sh was copied to "$simFoldDefaultFiles"
if [ -z $singleSiteName ]  # returns TRUE if unset or empty -> $siteList is used
 then # job array #
  array="T"   # needed by cleanup_slm_hpc_run.sh
  # get number of lines in the file (number of sites stored in the $siteList = $siteNumber -1 and starts at line 2, see jobArray below)
  siteNumber=$(wc -l "$siteList" | awk '{print $1}')
  # overwrite siteNumber with number of parameter combinations for parameter sensitivity
  if [ $paramSensitivitySim == "T" ]; then
    siteNumber=$(($paramSensi_NumberCombinations +1)) # jobArray generally starts with 2 and goes to number+1
  fi
  # compile job array argument 
  jobArray=$(echo "$jobName"[2-"$siteNumber"]%"$numberCpuMax")
  # submit job array starting at 2 because first line (in the all_site_list.dat) is the header
  bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobArray" -o "$simFoldOut"out.%I.txt -e "$simFoldErr"err.%I.txt "$emailSwitch" "$simFoldDefaultFiles""$slmRunScriptJobFilename" "$jobName" "$simDefaultFiles"
else # job with individual site #
  array="F"   # needed by cleanup_slm_hpc_run.sh
  bsub -q "$queue" -M "$memory" -R "$reservedArg" -J "$jobName" -o "$simFoldOut"out.txt -e "$simFoldErr"err.txt "$emailSwitch" "$simFoldDefaultFiles""$slmRunScriptJobFilename" "$jobName" "$simDefaultFiles" "$singleSiteName"
fi

##################################################################################
## cleanup temporary files and folders via cleanup_slm_hpc_run.sh               ##
##################################################################################
# this job is pending until the above started $jobName has ended
# the memory requirement is set to 1GB
sleep 1s  # just some waiting time to avoid any interference with the above started job
cleanupJobName="$jobName""Cleanup"
cleanQueue="short"
# start the job
bsub -w "ended($jobName)" -q "$cleanQueue" -M "1000000" -R "rusage[mem=1000]" -J "$cleanupJobName" -o "$simFold"cleanup.eo "$emailSwitch" "$simFoldDefaultFiles""$slmCleanupScriptFilename" "$jobName" "$array" "$singleSiteName"


## end prep_run_slm.sh ##
