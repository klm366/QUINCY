#!/bin/bash

#########
## Where: you can start this script from anywhere (no need to cd into a specific directory), when logged into a cluster node such as pc028
## What:  this script starts the SLM binary land.x , the specific model may be specified via the "slmModel" option
## HowTo: specify the below simulation options as suitable for you
## Tip:   some of the below simulation options can be provided as argument to this script, for more information type: ./run-slm_test-bed.sh -h 
## 
#########

########################################################################################################
####### namelist options you might want to change ######################################################
########################################################################################################

# model
slmModel="land"   # specific SLM model [default: "land"]; possible values: land plant soil canopy test_canopy test_radiation

# flags: T (True) or F (False) [large caps]
# plant-soil model options
sb_model_scheme="simple_0d"     # simple_0d simple_1d jsm 
sb_nloss_scheme="dynamic"       # none fixed dynamic 
flag_sb_prescribe_nh4="F"       # set nh4solute to a constant value
flag_sb_prescribe_no3="F"       # set no3solute to a constant value
flag_sb_prescribe_po4="T"       # set po4solute to a constant value
n_deposition_scheme="none"      # options: none constant dynamic (nhx_srf_down & noy_srf_down)
p_deposition_scheme="none"      # options: none constant dynamic (p_srf_down)

# BNF & Mycorrhiza
bnf_scheme="standard"           # none fixed standard optimal resistance unlimited
flag_mycorrhiza="F"             # TRUE: the other flag_myco* can be set TRUE | FALSE: the other flag_myco* are set to FALSE by this script
flag_mycorrhiza_org="F"         # allows mycorrhiza to access the organic/residual N pool
flag_mycorrhiza_prim="F"        # present mycorrhiza enhance decomposition of residual pool (priming)
# canopy model options
flag_t_resp_acclimation="T"
flag_t_jmax_acclimation="T"
flag_optimal_Nfraction="F"
ncanopy=10                      # number of canopy layers
canopy_layer_scheme="fapar"     # canopy layer depth calculation options: standard fapar

# plant growth model options
biomass_alloc_scheme="fixed"    # fixed dynamic optimal
leaf_stoichom_scheme="fixed"    # fixed dynamic optimal
# vegetation dynamics options
veg_dynamics_scheme="population"  # none population cohort
cohort_harvest_interval=80        # number of years after which the cohort is harvested and re-planted, is used with veg_dynamics_scheme="cohort"

# soil model options
soil_depth=5.7                  # meter (5.7 to 9.5)
nsoil_layer=5                   # number of soil layers (default for SLM land: 5 | for SLM soil only available forcing: 15 | generally: min 5 & max 20)

# JSM
flag_sb_jsm_transport="F"       # ..
flag_sb_jsm_litter_input="F"    # ..
flag_sb_jsm_OM_sorption="F"     # ..

# ECA and P cycle  (currently not available in the prep_run_slm.sh; defaults are: T/T/eca_full)
flag_sb_eca_upt_po4="T"         # TRUE: use the ECA approach for plant and microbe P uptake | FALSE: use non-ECA approach for plant and microbe P uptake
flag_sb_cle_po4="T"             # TRUE: Langmuir equilibrium is constantly maintained between solute P and labile P | FALSE: Langmuir equilibrium is not maintained
sb_adsorp_scheme="eca_full"     # eca_full; eca_part (only uptake, not for adsorption); eca_none

# simulation length & output
simulation_length_number=20     # 1 to ...500 ? [simulation_length_unit defaults to "y" == years]
output_interval_pool="weekly"   # timestep daily weekly monthly yearly
output_interval_flux="weekly"   # timestep daily weekly monthly yearly
output_start_first_day_year=1   # 1 to ... [lower than or equal with output_end_last_day_year ]
output_end_last_day_year=20     # equal with or larger than output_start_first_day_year

# PFT 
setPFTAutomatically="yes"       # yes / no (in case of 'yes', the value of option plant_functional_type_id is ignored)
plant_functional_type_id=5      # 1 - 8

# forcing input
flag_oneYearForcing="F"         # T (True) or F (False): whether to use (repeatedly) only the 1st year of the climate forcing input file
forcingDataSet="fluxnet"        # options: fluxnet cruncep face manip fert craine photacc none (for 'none' provide: test_bed/climate.dat [land plant canopy] or soil_forcing_litterfall.dat & soil_forcing_soilphysics.dat [soil])
siteName="DE-Hai"               # get the Site-ID from the 1st column of the particular "site list .dat"
                                # examples: [FLUX] "DE-Hai" "FI-Hyy" "ZA-Kru"  [CRUN] "BHS_01" "HAI_01" "ILO_05"  
                                #           [FACE] "DK-AMB" "OR-ELE"  [MANIP] "ES-MaG" "ES-MaT"  [FERT] "HAR-CON_4" "HAR-HIG_4" "HAR-LOW_4" "HAR-CON_5"
                                #           [CRAINE] "CRA_01" "CRA_02" "CRA_18"  [PHOTACC] "ARC_01" "MON_02" "STR_01"
soil_awc_prescribe=300.0        # float !
soil_theta_prescribe=1.0        # float !
longitude=0.0                   # float !
latitude=50.0                   # float !


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################



## change into the directory where this script is located; thus, you can start this script from anywhere ##
cd "$( dirname "$0" )"


## deal with options and arguments passed to the script ##
# getopts
while getopts ":p:m:i:f:s:y:b:z:h" opt; do
  case "$opt" in
  "p")
    echo "   plant_functional_type_id=""$OPTARG"
    echo "   setPFTAutomatically is set to: no"
    plant_functional_type_id="$OPTARG"
    setPFTAutomatically="no"
    ;;
  "m")
    echo "   slmModel=""$OPTARG"
    slmModel="$OPTARG"
    ;;
  "i")
    echo "   siteName=""$OPTARG"
    siteName="$OPTARG"
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
    "no")
      echo "   forcingDataSet=none"
      forcingDataSet="none"
      ;;
    esac
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
      echo "     mycorrhiza not active!"
      flag_mycorrhiza="$OPTARG"
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
    echo " -p set the PFT (plant_functional_type_id)"
    echo "    Arguments to -p:"
    echo "       any integer value 1-8"
    echo "       note: setPFTAutomatically is set to: no"
    echo " -m set slmModel"
    echo "    Arguments to -m:"
    echo "       land plant soil canopy test_canopy test_radiation"
    echo " -i set siteName, i.e., siteID"
    echo "    Arguments to -i:"
    echo "       a valid siteName, e.g., [FLUX] DE-Hai  [CRUN] BHS_01  [FACE] DK-AMB  [MANIP] ES-MaG  [FERT] HAR-CON_4  [CRAINE] CRA_01  [PHOTACC] MON_02"
    echo " -f set the forcingDataSet"
    echo "    Arguments to -f:"
    echo "       fl: FLUXNET"
    echo "       cr: CRUNCEPv7"
    echo "       fa: FACE sites"
    echo "       ma: MANIP tree and grassland sites"
    echo "       fe: FERT sites (fertilization)"
    echo "       ca: CRAINE sites (selected sites from doi:10.5061/dryad.v2k2607)"
    echo "       ph: PHOTACC sites (phot acclimation)"
    echo "       no: none; using climate.dat file in test_bed folder"
    echo " -s set the soil model"
    echo "    Arguments to -s:"
    echo "       s0d: simple_0d"
    echo "       s1d: simple_1d"
    echo "       jsm: jsm"    
    echo " -y set simulation length"
    echo "    Arguments to -y:"
    echo "       any integer value"
    echo "       note: this sets simulation output accordingly (1-x)"
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
    exit 1
    ;;
  esac
done

## 
echo " "
echo -n " Processing your request, please hold ... "

## prepare Validation of user-defined model-options ##
# Function to validate model-options provided by the user
# $1 must be a part of $2 list, else the script is stopped
# example: validateUserConfigValues "$queue" "$validateQueue"
function validateUserConfigValues {
validate="$(grep -w "$1" <<< "$2")"
if [ -z "${validate// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  echo "EXIT! > ""$1"" <  is not a valid value for parameter ""$3"", please change and restart run-slm_test-bed.sh"
  exit
 else
  if [ -z $1 ]
   then
    echo "EXIT!  The value for parameter ""$3"" is empty but must contain a valid string, please change and restart run-slm_test-bed.sh"
    exit
  fi
fi
return 0
}
# lists to validate against
validateSlmModel="land plant soil canopy test_canopy test_radiation"
validateBNFScheme="none fixed standard optimal resistance unlimited"
validateMycScheme="no am ecm"
validateNDepositionScheme="none constant dynamic"
validatePDepositionScheme="none constant dynamic"
validateCanopyLayerScheme="standard fapar"
validateSBModelScheme="simple_0d simple_1d jsm"
validateSBNLossScheme="none fixed dynamic"
validateBiomassAllocScheme="fixed dynamic optimal"
validateLeafStoichomScheme="fixed dynamic optimal"
validateVegDynamicsScheme="none population cohort"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validatePlantFunctionalTypeID="1 2 3 4 5 6 7 8"
validateOptionForcing="constant generate readfile"
validateSetPFTAutomatically="yes no"
validateForcingDataSet="fluxnet cruncep face manip fert craine photacc none"
validateSb_adsorp_scheme="eca_full eca_part eca_none"

## validate user-defined model-options ##
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$bnf_scheme" "$validateBNFScheme" "bnf_scheme"
validateUserConfigValues "$flag_mycorrhiza" "$validateFlags" "flag_mycorrhiza"
validateUserConfigValues "$flag_mycorrhiza_org" "$validateFlags" "flag_mycorrhiza_org"
validateUserConfigValues "$flag_mycorrhiza_prim" "$validateFlags" "flag_mycorrhiza_prim"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$sb_model_scheme" "$validateSBModelScheme" "sb_model_scheme"
validateUserConfigValues "$sb_nloss_scheme" "$validateSBNLossScheme" "sb_nloss_scheme"
validateUserConfigValues "$veg_dynamics_scheme" "$validateVegDynamicsScheme" "veg_dynamics_scheme"
validateUserConfigValues "$flag_sb_prescribe_nh4" "$validateFlags" "flag_sb_prescribe_nh4"
validateUserConfigValues "$flag_sb_prescribe_no3" "$validateFlags" "flag_sb_prescribe_no3"
validateUserConfigValues "$flag_sb_prescribe_po4" "$validateFlags" "flag_sb_prescribe_po4"
validateUserConfigValues "$flag_t_resp_acclimation" "$validateFlags" "flag_t_resp_acclimation"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_sb_jsm_transport" "$validateFlags" "flag_sb_jsm_transport"
validateUserConfigValues "$flag_sb_jsm_litter_input" "$validateFlags" "flag_sb_jsm_litter_input"
validateUserConfigValues "$flag_sb_jsm_OM_sorption" "$validateFlags" "flag_sb_jsm_OM_sorption"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$biomass_alloc_scheme" "$validateBiomassAllocScheme" "biomass_alloc_scheme"
validateUserConfigValues "$leaf_stoichom_scheme" "$validateLeafStoichomScheme" "leaf_stoichom_scheme"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$setPFTAutomatically" "$validateSetPFTAutomatically" "setPFTAutomatically"
validateUserConfigValues "$plant_functional_type_id" "$validatePlantFunctionalTypeID" "plant_functional_type_id"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$forcingDataSet" "$validateForcingDataSet" "forcingDataSet"
validateUserConfigValues "$flag_sb_eca_upt_po4" "$validateFlags" "flag_sb_eca_upt_po4"
validateUserConfigValues "$flag_sb_cle_po4" "$validateFlags" "flag_sb_cle_po4"
validateUserConfigValues "$sb_adsorp_scheme" "$validateSb_adsorp_scheme" "sb_adsorp_scheme"

# not validated: simulation_length_number, output_start_first_day_year, output_end_last_day_year, soil_awc_prescribe, soil_theta_prescribe, latitude, longitude 


## definition of some basic variables
siteDataDirBase="/Net/Groups/BGI/data/Model/SLM1/input/point/"
metForcingDir="met_forcing/"
soilForcingDirBase="soil_forcing_sl"
forceSoilLitterFile="soil_forcing_litterfall.dat"
forceSoilPhysicFile="soil_forcing_soilphysics.dat"
forceLandMetFile="climate.dat"


## set folder for the site data-set
case "$forcingDataSet" in
"fluxnet")
  siteDataSet="FLUXNET2/"
  ;;
"cruncep")
  siteDataSet="CRUNCEPv7/"
  ;;
"face")
  siteDataSet="FACE/"
  ;;
"manip")
  siteDataSet="MANIP/"
  ;;
"fert")
  siteDataSet="FERT/"
  ;;
"craine")
  siteDataSet="CRAINE/"
  ;;
"photacc")
  siteDataSet="PHOTACC/"
  ;;
esac
# no soil forcing for MANIP FERT CRAINE PHOTACC
if [ "$slmModel" == "soil" ]; then
  if [ "$forcingDataSet" == "manip" ] || [ "$forcingDataSet" == "fert" ] || [ "$forcingDataSet" == "craine" ] || [ "$forcingDataSet" == "photacc" ]; then
    echo " Sorry, no soil forcing available for ""$forcingDataSet"" sites yet."
    echo " EXIT!"
    exit 1
  fi
fi


## test number of soil levels and set dirname for soil_forcing files
if [ $slmModel == "soil" ]
  then
  # test number of soil layers against the default value 15; soil-forcing is available for 15 sl only
  if [ "$nsoil_layer" != 15 ]; then
    echo " "
    echo " You aim at using ""$nsoil_layer"" soil layers, where soil forcing is available for 15 sl only."
    echo " nsoil_layer is now set to 15."
    nsoil_layer=15
  fi
  # test soil depth against the value necessary to run 15 sl, i.e., 9.5 meters
  if [ "$soil_depth" != 9.5 ]; then
    echo " You aim at using 15 soil layers, which needs a soil depth of 9.5 meters. But your soil depth is set to ""$soil_depth""."
    echo " soil_depth is now set to 9.5 meters."
    soil_depth=9.5
  fi
    # soil forcing directory
    soilForcingDir="$soilForcingDirBase""15/"
fi 


## define site info file (full path)
case "$forcingDataSet" in
"fluxnet")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"flux_all_sites_1995-2006_list.dat
  ;;
"cruncep")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"crun_cue_gfdb_sites_1984-2013_list.dat
  # the two cruncep spp sites are not listed in the all_gfdb_sites_*_list.dat but only in the spp_sites_*_list.dat
  if [ "$siteName" == "BBR_01" ] || [ "$siteName" == "LUE_01" ] || [ "$siteName" == "VES_01" ]
   then
    allSiteInfoFile="$siteDataDirBase""$siteDataSet"crun_spp_sites_1984-2013_list.dat
  fi
  ;;
"face")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"face_sites_1996-2008_list.dat
  ;;
"manip")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"manip_sites_2009-2017_list.dat
  ;;
"fert")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"fert_sites_1901-2015_list.dat
  echo " "
  echo " Note: you are running (a) FERT site(s) using forcing 1901-2015 by default."
  ;;
"craine")
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"craine_sites_1901-2015.dat
  echo " "
  echo " Note: you are running (a) CRAINE site(s) using forcing 1901-2015 by default."
  ;;
"photacc")
    #! temporary fix for PHOTACC !#
    siteDataDirBase="/Net/Groups/BGI/scratch/jengel/quincy_data/point/"
  allSiteInfoFile="$siteDataDirBase""$siteDataSet"photacc_sites_1901-2015.dat
  echo " "
  echo " Note: you are running (a) PHOTACC site(s) using forcing 1901-2015 by default."
  ;;
esac


## test if siteName exists in site info file
if [ "$forcingDataSet" != "none" ]; then
  # validate if $siteName variable has some (any) content
  if [ -z "$siteName" ]; then
    echo -e '\n'"  You did not set a site name. Variable siteName is empty."
    echo "EXIT!"
    exit 1
  fi
  # perform search for siteName in site list
  allSites=$(cat "$allSiteInfoFile" | awk '{print $1}')
  validateUserConfigValues "$siteName" "$allSites" "siteName"
  # in addition: show a note for cruncep site limitation of soil forcing available to cue and spp sites, but do not explicitly test for that
  if [ $slmModel == "soil" ] && [ "$forcingDataSet" == "cruncep" ]; then
    echo -e '\n'"  Your are using SLM soil for a cruncep site."
    echo    "  Note, soil forcing is available only for cue and spp sites."
  fi
fi


## get an identifier for the met_forcing file name for static forcing
#  i.e. the met_forcing file for the specific site (not the full file name but an string long enough to identify the file correctly)
case "$forcingDataSet" in
"fluxnet")
  siteNameStatMetForc="$siteName""_s_"                # e.g.: DE-Hai  ->  DE-Hai_s_*.dat  ->  DE-Hai_s_2000-2006.dat  -> climate.dat
  ;;
"cruncep")
  crunSiteNameShort=$(echo "$siteName" | cut -c1-3)   # shortens file name provided, e.g.: BHS_01 ->  BHS  (cruncep met_forcing is not sub-site specific)
  siteNameStatMetForc="$crunSiteNameShort""_1984"
  ;;
"face")
  siteNameStatMetForc="$siteName""_199"               # selects the static focing file such as DK-AMB_1996-2007 but not the transient DK-AMB_1901-2015
  ;;
"manip")
  siteNameStatMetForc="$siteName""_2009"              # ES-MaG_2009-2017.dat & ES-MaT_2009-2017.dat
  ;;
"fert")
  siteNameStatMetForc="$siteName""_1901"              # HAR-CON_1901-2015.dat (only transient forcing files)
  ;;
"craine")
  siteNameStatMetForc="$siteName""_1901"              # CRA_01_1901-2015.dat (only transient forcing files)
  ;;
"photacc")
  siteNameStatMetForc="$siteName""_1901"              # DRE_01_1901-2015.dat (only transient forcing files)
  ;;
esac


## test_canopy/test_radiation model and forcingDataSet none ##
if [ $slmModel != "test_canopy" ] && [ $slmModel != "test_radiation" ]; then
  # no site data-set
  if [ "$forcingDataSet" == "none" ]  # if forcing data are neither Fluxnet nor Cruncep dataset
  then  # local forcing file
    # SLM soil / land plant canopy
    if [ $slmModel == "soil" ]
      then
        echo -e '\n' "  Info: you set forcingDataSet to 'none' using the test_bed/soil_forcing_litterfall.dat & soil_forcing_soilphysics.dat files."
        echo "   No site specific parameter values are available in this mode, and the option setPFTAutomatically is used with 'no'."
        if [ ! -e ./"$forceSoilLitterFile" ] # if file not exists
          then
            echo "EXIT!  the test_bed/""$forceSoilLitterFile"" does not exist, but the option forcingDataSet is set to none !"
            exit
        fi
        if [ ! -e ./"$forceSoilPhysicFile" ] # if file not exists
          then
            echo "EXIT!  the test_bed/""$forceSoilPhysicFile"" does not exist, but the option forcingDataSet is set to none !"
            exit
        fi
    
    else # land plant canopy
      echo -e '\n' "  Info: you set forcingDataSet to 'none' using the test_bed/climate.dat."
      echo "   No site specific parameter values are available in this mode, and the option setPFTAutomatically is used with 'no'."
      echo "   PFT is set to: ""$plant_functional_type_id"
      if [ ! -e ./"$forceLandMetFile" ] # if file not exists
        then
          echo "EXIT!  the test_bed/""$forceLandMetFile"" does not exist, but the option forcingDataSet is set to none !"
          exit
      fi
    fi 
    setPFTAutomatically="no"   # for forcingDataSet = none
  fi  # end: forcingDataSet == "none"
else  
  setPFTAutomatically="no"   # for test_canopy/test_radiation model
fi


## create folder, and copy some files there ##
# postprocessing R scripts
# namelist
# lctlib_slm_nlct8.def
mkdir -p ./"$slmModel"
rm -rf   ./"$slmModel"/*  # delete all in this folder if the folder existed before
cp ../r_script/postprocessing_"$slmModel".R   ./"$slmModel"
cp ../namelist/namelist.slm                   ./"$slmModel"
cp ../../data/lctlib_slm_nlct8.def            ./"$slmModel"


## copy forcing file(s)
if [ $slmModel != "test_canopy" ] && [ $slmModel != "test_radiation" ]; then      # copy a climate.dat file for all models but test_canopy & test_radiation
    # forcingDataSet == none (local forcing file needed)
    if [ "$forcingDataSet" == "none" ]
      then
        if [ $slmModel == "soil" ]
          then
            cp ./"$forceSoilLitterFile" ./"$slmModel"/"$forceSoilLitterFile"
            cp ./"$forceSoilPhysicFile" ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp ./"$forceLandMetFile"    ./"$slmModel"/"$forceLandMetFile"
        fi 
      # forcingDataSet == fluxnet / cruncep / face / manip / fert / craine / photacc
      else 
        if [ $slmModel == "soil" ]
          then
            # use a shorter site name for cruncep forcing 
            if [ "$forcingDataSet" == "cruncep" ]; then
              siteNameSoilForceID=$(echo "$siteName" | cut -c1-3)   # e.g. ILO_05 -> ILO
            else # fluxnet / face / manip / fert / craine / photacc
              siteNameSoilForceID="$siteName"
            fi
            # copy files
            cp "$siteDataDirBase""$siteDataSet""$soilForcingDir""$siteNameSoilForceID"*_lit.dat  ./"$slmModel"/"$forceSoilLitterFile"
            cp "$siteDataDirBase""$siteDataSet""$soilForcingDir""$siteNameSoilForceID"*_sph.dat  ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp "$siteDataDirBase""$siteDataSet""$metForcingDir""$siteNameStatMetForc"*.dat  ./"$slmModel"/"$forceLandMetFile"
        fi 
    fi
fi


## test cohort-harvesting interval, if veg_dynamics_scheme="cohort"
if [ "$veg_dynamics_scheme" == "cohort" ]
 then
   if [ "$cohort_harvest_interval" -lt 20 ] || [ "$cohort_harvest_interval" -gt 140 ]
    then
      echo " Note. The cohort-harvesting interval might be set to a value from 20 to 140. You set it to: ""$cohort_harvest_interval""."
   fi
   if [ "$cohort_harvest_interval" -lt 1 ]
    then
      echo " The cohort-harvesting interval is set to a value smaller than 1."
      echo " EXIT!"
      exit 1
   fi
fi


## set mycorrhiza flag_myco_* flags to FALSE if mycorrhiza are disabled overall
if [ $flag_mycorrhiza == "F" ]
 then
  flag_mycorrhiza_org="F"
  flag_mycorrhiza_prim="F"
fi


## cd into the new folder ##
cd ./"$slmModel"

## manipulate namelist.slm accordingly to above options ##
 # function: pass the variable_name ($1) & the new_value ($2)
 # this function assumes a specific relative path to the namlist.slm to be valid
function changeNamelistValue {
oldValue=$(grep -w $1 ./namelist.slm | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 ./namelist.slm | cut -d : -f 1)                       # get line number of variable $1
if [[ $(uname -s) == "Darwin" ]]; then                                      # make it work at Notebook "Darwin" -> Mac operating system
  gsed -i "$lineNum"s/$oldValue/$2/ ./namelist.slm                          # replace the old value in line $lineNum with $2
else 
  sed -i "$lineNum"s/$oldValue/$2/ ./namelist.slm
fi
#sync
return 0                                                                              # get return value with: return_value=$?
}
 # examples:
 #changeNamelistValue simulation_length_number 300
 #changeNamelistValue output_interval_flux monthly    # quotation marks not needed
 #changeNamelistValue flag_calc_soil_moisture .FALSE.
 #changeNamelistValue soil_theta_prescribe 1.

changeNamelistValue slm_model_name "$slmModel"
changeNamelistValue bnf_scheme "$bnf_scheme" 
changeNamelistValue n_deposition_scheme "$n_deposition_scheme"
changeNamelistValue p_deposition_scheme "$p_deposition_scheme"
changeNamelistValue ncanopy "$ncanopy"
changeNamelistValue canopy_layer_scheme "$canopy_layer_scheme"
changeNamelistValue veg_dynamics_scheme "$veg_dynamics_scheme"
changeNamelistValue biomass_alloc_scheme "$biomass_alloc_scheme"
changeNamelistValue leaf_stoichom_scheme "$leaf_stoichom_scheme"
changeNamelistValue cohort_harvest_interval "$cohort_harvest_interval"
changeNamelistValue sb_model_scheme "$sb_model_scheme"
changeNamelistValue sb_nloss_scheme "$sb_nloss_scheme"
changeNamelistValue sb_bnf_scheme "$bnf_scheme" 
changeNamelistValue sb_adsorp_scheme "$sb_adsorp_scheme"
if [ $flag_sb_eca_upt_po4 == T ]
  then
    changeNamelistValue flag_sb_ECA_upt_po4 ".TRUE."
  else
    changeNamelistValue flag_sb_ECA_upt_po4 ".FALSE."
fi
if [ $flag_sb_cle_po4 == T ]
  then
    changeNamelistValue flag_sb_CLE_po4 ".TRUE."
  else
    changeNamelistValue flag_sb_CLE_po4 ".FALSE."
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
if [ $flag_oneYearForcing == T ]
  then
    changeNamelistValue flag_climate_forcing_one_year ".TRUE."
  else
    changeNamelistValue flag_climate_forcing_one_year ".FALSE."
fi
changeNamelistValue simulation_length_number "$simulation_length_number"
changeNamelistValue output_interval_pool "$output_interval_pool"
changeNamelistValue output_interval_flux "$output_interval_flux"
changeNamelistValue output_start_first_day_year "$output_start_first_day_year"
changeNamelistValue output_end_last_day_year "$output_end_last_day_year"
changeNamelistValue soil_awc_prescribe "$soil_awc_prescribe"
changeNamelistValue soil_theta_prescribe "$soil_theta_prescribe"
changeNamelistValue soil_depth "$soil_depth"
changeNamelistValue nsoil_energy "$nsoil_layer"
changeNamelistValue nsoil_water "$nsoil_layer"
changeNamelistValue set_parameter_values_from_file ".FALSE."   # set to default FALSE because parameter update is not applicable here in test_bed

## set site-specific information (not for test_canopy nor test_radiation)
if [ "$forcingDataSet" != "none" ] && [ "$slmModel" != "test_canopy" ] && [ "$slmModel" != "test_radiation" ]
  then
    # get data of the particular site
    siteData=$(grep "$siteName" "$allSiteInfoFile")
    # change namelist values accordingly, column numbers from all_sites_list.dat are:
    # 1:Site-ID 2:lon 3:lat 4:PFT 5:start 6:end 7:lon_gmt 8:clay 9:silt 10:sand 
    # 11:awc 12:bd 13:pH 14:taxusda 15:taxnwrb 16:lith_glim 17:LAI 18:Nleaf 19:SLA 20:Height 
    # 21:PlantYear 22:soil_p_depth 23:soil_p_labile 24:soil_p_slow 25:soil_p_occluded 26:soil_p_primary 27:qmax_org_fine_particle 
    changeNamelistValue longitude               "$(awk '{print $2}'  <<< "$siteData")"
    changeNamelistValue latitude                "$(awk '{print $3}'  <<< "$siteData")"
    pftName="$(awk '{print $4}' <<< "$siteData")"
    changeNamelistValue soil_clay               "$(awk '{print $8}'  <<< "$siteData")"
    changeNamelistValue soil_silt               "$(awk '{print $9}'  <<< "$siteData")"
    changeNamelistValue soil_sand               "$(awk '{print $10}' <<< "$siteData")"
    changeNamelistValue bulk_density            "$(awk '{print $12}' <<< "$siteData")"
    changeNamelistValue soil_ph                 "$(awk '{print $13}' <<< "$siteData")"
    changeNamelistValue usda_taxonomy_class     "$(awk '{print $14}' <<< "$siteData")"
    changeNamelistValue nwrb_taxonomy_class     "$(awk '{print $15}' <<< "$siteData")"
    changeNamelistValue lithology_class         "$(awk '{print $16}' <<< "$siteData")"
    changeNamelistValue lai_site_specific       "$(awk '{print $17}' <<< "$siteData")"
    changeNamelistValue soil_p_depth            "$(awk '{print $22}' <<< "$siteData")"
    changeNamelistValue soil_p_labile           "$(awk '{print $23}' <<< "$siteData")"
    changeNamelistValue soil_p_slow             "$(awk '{print $24}' <<< "$siteData")"
    changeNamelistValue soil_p_occluded         "$(awk '{print $25}' <<< "$siteData")"
    changeNamelistValue soil_p_primary          "$(awk '{print $26}' <<< "$siteData")"
    changeNamelistValue qmax_org_fine_particle  "$(awk '{print $27}' <<< "$siteData")"
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
        echo "EXIT!  no such PFT is known (failure during trying to automatically grep the PFT ID)."
        exit 
    esac
    # set PFT
    if [ "$setPFTAutomatically" == "yes" ]
     then
      changeNamelistValue plant_functional_type_id "$pftID"
    else   
      changeNamelistValue plant_functional_type_id "$plant_functional_type_id"
    fi

  else  # forcingDataSet == none || test_canopy || test_radiation
    changeNamelistValue plant_functional_type_id "$plant_functional_type_id"
    changeNamelistValue longitude "$longitude"
    changeNamelistValue latitude "$latitude"
fi


## confirmation of finishing this script & starting SLM ##
echo " Nice idea !" 
if [ "$forcingDataSet" == "none" ]; then
  # siteName variable not valid, so give different feedback compared to using a forcing data set
  echo " Let's go with SLM ""$slmModel"" for an unknown site, using the local forcing file." 
else
  echo " Let's go with SLM ""$slmModel"" for the site: ""$siteName" 
fi
echo -n "  ...  "

## run the SLM binary (always compiled as land.x) ##
compiler=$(sed -n 2p ./../../../build/Make.inc.local | awk '{print $3}')
./../../../"$compiler"/bin/land.x 
echo "Simulation finished"
echo " "
echo "Running postprocessing ... "
#sync

## run R script for postprocessing ##
Rscript --vanilla ./postprocessing_"$slmModel".R

#sync
echo " "
echo " ... done !"

## end of this script, confirmation about successful simulation is provided by the R script itself ##
