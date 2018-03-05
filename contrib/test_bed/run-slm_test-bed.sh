#!/bin/bash

#########
## Where: you can start this script from anywhere (no need to cd into a specific directory), when logged into a cluster node such as pc028
## What:  this script starts the SLM binary land.x , the specific model may be specified via the "slmModel" option
## HowTo: specify the below simulation options as suitable for you
## Tip:   some of the below simulation options can be provided as argument to this script, for more information type: ./run-slm_test-bed.sh -h 
## 
## Note:  the namelist option 'set_parameter_values_from_file' is set to .FALSE. by this script - please use the pre_run.sh @ HPC (Cluster) if you need .TRUE.

########################################################################################################
####### namelist options you might want to change ######################################################
########################################################################################################

# model
slmModel="land"   # specific SLM model [default: "land"]; possible values: land plant soil canopy test_canopy

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
flag_myco_inorganic="F"         # allows mycorrhiza to access the inorganic N pools
flag_myco_organic="F"           # allows mycorrhiza to access the organic/residual N pool
flag_myco_priming="F"           # present mycorrhiza enhance decomposition of residual pool (priming)
# canopy model options
flag_t_resp_acclimation="T"
flag_t_jmax_acclimation="T"
flag_optimal_Nfraction="F"
ncanopy=10                      # number of canopy layers
canopy_layer_scheme="fapar"     # canopy layer depth calculation options: standard fapar

# plant growth model options
flag_optimal_Nleaf="F"
flag_optimal_allocation="F"
# vegetation dynamics options
flag_tile_veg_dynamics="T"      # does also set "start_bareground"

# soil model options
soil_depth=9.5                  # meter
nsoil_layer=15                   # number of soil layers (default: 5 | min: 5 | max: 20)

# ECA and P cycle
flag_sb_eca_upt_po4="T"         # TRUE: use the ECA approach for plant and microbe P uptake | FALSE: use non-ECA approach for plant and microbe P uptake
flag_sb_cle_po4="T"             # TRUE: Langmuir equilibrium is constantly maintained between solute P and labile P | FALSE: Langmuir equilibrium is not maintained
flag_sb_lk_po4="T"              # TRUE: Langmuir kinetics is used for adsorption flux | FALSE: adsorption rate is determined by CLE
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
forcingDataSet="fluxnet"        # options: fluxnet cruncep face none (for 'none' provide: test_bed/climate.dat [land plant canopy] or soil_forcing_litterfall.dat & soil_forcing_soilphysics.dat [soil])
siteName="DE-Hai_2000-2006"     # provide filename (without the trailing '.dat') of any site from /Net/Groups/BGI/data/Model/SLM1/input/point/.. FLUXNET2/ or CRUNCEPv7/ or FACE/ 
                                # e.g. [FLUX] DE-Hai_2000-2006, AT-Neu_2002-2005 [CRUN] HOW_1984-2013, KAN_1984-2013 [FACE] DK-AMB_1996-2007, OR-AMB_1998-2008
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
# check if '-o' option is used, if so, set all optimality flags to F per default
argOpt_o="$(grep -w "\-o" <<< "$@")"
if [ ! -z "${argOpt_o// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  flag_optimal_Nfraction="F"
  flag_optimal_Nleaf="F"
  flag_optimal_allocation="F"
fi
# getopts
while getopts ":p:m:s:y:o:b:z:h" opt; do
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
  "o")
    case "$OPTARG" in
    "no")
      echo "   no optimality"
      flag_optimal_Nfraction="F"
      flag_optimal_Nleaf="F"
      flag_optimal_allocation="F"
      ;;
    "all")
      echo "   using optimal Nfraction"
      echo "   using optimal Nleaf"      
      echo "   using optimal Allocation"
      flag_optimal_Nfraction="T"
      flag_optimal_Nleaf="T"
      flag_optimal_allocation="T"
      ;;
    "nf")
      echo "   using optimal Nfraction"
      flag_optimal_Nfraction="T"
      ;;
    "nl")
      echo "   using optimal Nleaf"
      flag_optimal_Nleaf="T"
      ;;
    "al")
      echo "   using optimal Allocation"
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
    echo " -o set the optimality scheme(s)"
    echo "    Arguments to -o:"
    echo "       no : no optimality"
    echo "       all: all optimality schemes"
    echo "       nf : using only optimal Nfraction"
    echo "       nl : using only optimal Nleaf"
    echo "       al : using only optimal allocation"
    echo "    note: you can use the -o option multiple times as argument to the script, e.g., script.sh -o nf -o al"
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
validateSlmModel="land plant soil canopy test_canopy"
validateBNFScheme="none fixed standard optimal resistance unlimited"
validateMycScheme="no am ecm"
validateNDepositionScheme="none constant dynamic"
validatePDepositionScheme="none constant dynamic"
validateCanopyLayerScheme="standard fapar"
validateSBModelScheme="simple_0d simple_1d jsm"
validateSBNLossScheme="none fixed dynamic"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validatePlantFunctionalTypeID="1 2 3 4 5 6 7 8"
validateOptionForcing="constant generate readfile"
validateSetPFTAutomatically="yes no"
validateForcingDataSet="fluxnet cruncep face none"
validateSb_adsorp_scheme="eca_full eca_part eca_none"

## validate user-defined model-options ##
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$bnf_scheme" "$validateBNFScheme" "bnf_scheme"
validateUserConfigValues "$flag_mycorrhiza" "$validateFlags" "flag_mycorrhiza"
validateUserConfigValues "$flag_myco_inorganic" "$validateFlags" "flag_myco_inorganic"
validateUserConfigValues "$flag_myco_organic" "$validateFlags" "flag_myco_organic"
validateUserConfigValues "$flag_myco_priming" "$validateFlags" "flag_myco_priming"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$n_deposition_scheme" "$validateNDepositionScheme" "n_deposition_scheme"
validateUserConfigValues "$p_deposition_scheme" "$validatePDepositionScheme" "p_deposition_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$sb_model_scheme" "$validateSBModelScheme" "sb_model_scheme"
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
validateUserConfigValues "$setPFTAutomatically" "$validateSetPFTAutomatically" "setPFTAutomatically"
validateUserConfigValues "$plant_functional_type_id" "$validatePlantFunctionalTypeID" "plant_functional_type_id"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$forcingDataSet" "$validateForcingDataSet" "forcingDataSet"
validateUserConfigValues "$flag_sb_eca_upt_po4" "$validateFlags" "flag_sb_eca_upt_po4"
validateUserConfigValues "$flag_sb_cle_po4" "$validateFlags" "flag_sb_cle_po4"
validateUserConfigValues "$flag_sb_lk_po4" "$validateFlags" "flag_sb_lk_po4"
validateUserConfigValues "$sb_adsorp_scheme" "$validateSb_adsorp_scheme" "sb_adsorp_scheme"

# not validated: simulation_length_number, output_start_first_day_year, output_end_last_day_year, soil_awc_prescribe, soil_theta_prescribe, latitude, longitude 


## definition of some basic variables
# model forcing - input files for the model
forceSoilLitterFile="soil_forcing_litterfall.dat"
forceSoilPhysicFile="soil_forcing_soilphysics.dat"
forceLandMetFile="climate.dat"

## test if forcing data exist in testbed/ or in the fluxnet/cruncep site-list-folder [if model is not test_canopy] ##
if [ $slmModel != "test_canopy" ] # copy a climate.dat file with PFT "BNE" for all models but test_canopy
  then
    if [ "$forcingDataSet" == "none" ]  # if forcing data are neither Fluxnet nor Cruncep dataset
    then  # local forcing file
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
      setPFTAutomatically="no"
    else  # fluxnet/cruncep forcing file
      case "$forcingDataSet" in
      "fluxnet")
        siteDataSetFull="FLUXNET2"
        ;;
      "cruncep")
        siteDataSetFull="CRUNCEPv7"
        ;;
      "face")
        siteDataSetFull="FACE"
        if [ $slmModel == "soil" ]
          then
            echo " Sorry, no soil_forcing files for FACE sites yet."
            echo " EXIT!"
            exit 1
        fi
        ;;
      esac

      # soil / met forcing
      if [ $slmModel == "soil" ]
        then
          #forcing="soil_forcing"
          if [ $nsoil_layer -eq 5 ]
            then
              forcing="soil_forcing_sl5"
          elif [ $nsoil_layer -eq 9 ]
            then
              forcing="soil_forcing_sl9"
          elif [ $nsoil_layer -eq 15 ]
            then
              forcing="soil_forcing_sl15"
          else
            echo "Note: soil forcing is available only for 5, 9 or 15 soil layers. You selected sl: ""$nsoil_layer""."
            echo "EXIT!"
            exit 1
          fi
        else # land plant canopy
          forcing="met_forcing"
      fi 

      siteFolderBasic="/Net/Groups/BGI/data/Model/SLM1/input/point/"
      siteFolder="$siteFolderBasic""$siteDataSetFull"
      siteFolderForcing="$siteFolder""/""$forcing"
      if ! ls -la $siteFolderForcing | grep $siteName > /dev/null  # if forcing file not found
        then
          echo "EXIT! The site ""$siteName"" you selected for forcing input does not exist in the forcingDataSet folder !"
          exit
      fi
  fi  # ["$forcingDataSet" == "none"]
else  # in case: $slmModel == "test_canopy"
  setPFTAutomatically="no"
fi    # [$slmModel != "test_canopy"]


## create folder, and copy some files there ##
# postprocessing R scripts
# namelist
# lctlib_slm_nlct8.def
# forcing (climate/ soil)
mkdir -p ./"$slmModel"
rm -rf   ./"$slmModel"/*  # delete all in this folder if the folder existed before
cp ../r_script/postprocessing_"$slmModel".R   ./"$slmModel"
cp ../namelist/namelist.slm                   ./"$slmModel"
cp ../../data/lctlib_slm_nlct8.def            ./"$slmModel"
if [ $slmModel != "test_canopy" ] # copy a climate.dat file for all models but test_canopy
  then
    if [ "$forcingDataSet" == "none" ]
      then
        if [ $slmModel == "soil" ]
          then
            cp ./"$forceSoilLitterFile" ./"$slmModel"/"$forceSoilLitterFile"
            cp ./"$forceSoilPhysicFile" ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp ./"$forceLandMetFile"    ./"$slmModel"/"$forceLandMetFile"
        fi 
      else # $forcingDataSet == fluxnet / cruncep / face
        if [ $slmModel == "soil" ]
          then
            cp "$siteFolderForcing"/"$siteName"_lit.dat ./"$slmModel"/"$forceSoilLitterFile"
            cp "$siteFolderForcing"/"$siteName"_sph.dat ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp "$siteFolderForcing"/"$siteName".dat     ./"$slmModel"/"$forceLandMetFile"
        fi 
    fi
fi


# test number of soil layers against min and max values (5:20)
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


## set mycorrhiza flag_myco_* flags to FALSE if mycorrhiza are disabled overall
if [ $flag_mycorrhiza == "F" ]
 then
  flag_myco_inorganic="F"
  flag_myco_organic="F"
  flag_myco_priming="F"
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
  gsed -i "$lineNum"s/$oldValue/$2/ ./namelist.slm                             # replace the old value in line $lineNum with $2
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
if [ $flag_sb_lk_po4 == T ]
  then
    changeNamelistValue flag_sb_LK_po4 ".TRUE."
  else
    changeNamelistValue flag_sb_LK_po4 ".FALSE."
fi
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

# forcingDataSet (not for test_canopy)
if [ "$forcingDataSet" != "none" ] && [ "$slmModel" != "test_canopy" ]
  then
    case "$forcingDataSet" in
    "fluxnet")
      allSiteDatFile="$siteFolder"/all_sites_list.dat
      siteNameShort=$(echo "$siteName" | cut -c1-6)             # shortens file name provided, e.g.: CA-Man_1999-2003 -> CA-Man
      ;;
    "cruncep")
      allSiteDatFile="$siteFolder"/all_gfdb_sites_1984-2013_list.dat
      siteNameShort=$(echo "$siteName" | cut -c1-3)             # shortens file name provided, e.g.: KAD_1901-2015 -> KAD
      siteNameShort="$siteNameShort""_01"                       # KAD -> KAD_01
      # the two cruncep spp sites are not listed in the all_gfdb_sites_*_list.dat but only in the spp_sites_*_list.dat
      if [ "$siteNameShort" == "BBR_01" ] || [ "$siteNameShort" == "LUE_01" ]
       then
        allSiteDatFile="$siteFolder"/spp_sites_1984-2013_list.dat
      fi
      ;;
    "face")
      allSiteDatFile="$siteFolder"/all_sites_1996-2008_list.dat
      siteNameShort=$(echo "$siteName" | cut -c1-6)             # shortens file name provided, e.g.: DK-AMB_1996-2007 -> DK-AMB
      ;;
    esac

    # get data of the particular site
    siteData=$(grep "$siteNameShort" "$allSiteDatFile")
    # change namelist values accordingly
    changeNamelistValue longitude             "$(awk '{print $2}'  <<< "$siteData")"
    changeNamelistValue latitude              "$(awk '{print $3}'  <<< "$siteData")"
    pftName="$(awk '{print $4}' <<< "$siteData")"
    changeNamelistValue soil_clay             "$(awk '{print $8}'  <<< "$siteData")"
    changeNamelistValue soil_silt             "$(awk '{print $9}'  <<< "$siteData")"
    changeNamelistValue soil_sand             "$(awk '{print $10}' <<< "$siteData")"
    changeNamelistValue bulk_density          "$(awk '{print $12}' <<< "$siteData")"
    changeNamelistValue soil_ph               "$(awk '{print $13}' <<< "$siteData")"
    changeNamelistValue usda_taxonomy_class   "$(awk '{print $14}' <<< "$siteData")"
    changeNamelistValue nwrb_taxonomy_class   "$(awk '{print $15}' <<< "$siteData")"
    changeNamelistValue lithology_class       "$(awk '{print $16}' <<< "$siteData")"
    changeNamelistValue lai                   "$(awk '{print $17}' <<< "$siteData")"
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

  else  # forcingDataSet == none
    changeNamelistValue plant_functional_type_id "$plant_functional_type_id"
    changeNamelistValue longitude "$longitude"
    changeNamelistValue latitude "$latitude"
fi


## confirmation of finishing this script & starting SLM ##
echo " Nice idea !" #OK"
echo " Let's go with SLM ""$slmModel" #Starting SLM Simulation: ""$slmModel"
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
