#!/bin/bash

## 20161212
## Where: you can start this script from anywhere (no need to cd into a specific directory), when logged into a cluster node such as pc028
## What:  this script starts the SLM binary land.x , the specific model may be specified via the "slmModel" option
## HowTo: specify the below options as suitable for you
## Tip:   you can pass the model name as argument to this script [i.e. run-slm_test-bed.sh modelName], which overwrites the below defined value of slmModel
## 
## Note:  the namelist option set_parameter_values_from_file is set to .FALSE. by this script - please use the pre_run.sh @ HPC (Cluster) if you need .TRUE.

########################################################################################################
####### namelist options you might want to change ######################################################
########################################################################################################

# model
slmModel="land"   # specific SLM model [default: "land"]; possible values: land plant soil canopy test_canopy

# flags: T (True) or F (False) [large caps]
# plant-soil model options
n_uptake_scheme="standard"      # prescribed standard optimal resistance
p_uptake_scheme="prescribed"    # prescribed standard 
soil_moisture_scheme="bucket"   # prescribed bucket layer
flag_calc_soil_energy="T"
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
soil_depth=5.7                  # meter

# simulation length & output
simulation_length_number=10     # 1 to ...500 ? [simulation_length_unit defaults to "y" == years]
output_interval_pool="weekly"   # timestep daily weekly monthly yearly
output_interval_flux="weekly"   # timestep daily weekly monthly yearly
output_start_first_day_year=1   # 1 to ... [lower than or equal with output_end_last_day_year ]
output_end_last_day_year=10     # equal with or larger than output_start_first_day_year

# PFT 
setPFTAutomatically="yes"       # yes / no (in case of 'yes', the value of option plant_functional_type_id is ignored)
plant_functional_type_id=5      # 1 - 8

# forcing input
flag_oneYearForcing="F"         # T (True) or F (False): whether to use (repeatedly) only the 1st year of the climate forcing input file
siteDataSet="f"                 # [f] FLUXNET [c] CRUNCEPv5 [n] none (in case of 'none' a test_bed/climate.dat is expected)
siteName="DE-Hai_2000-2006"     # filename of any site from /Net/Groups/BGI/data/Model/SLM1/input/point/.. FLUXNET/ or CRUNCEPv5/ without the trailing '.dat'
                                # e.g. FLUX: DE-Hai_2000-2006, AT-Neu_2002-2005, ES-ES1_1999-2004 & CRUN: THO_D15_1984-2013, CAL_1984-2013, UNI_1984-2013
soil_awc_prescribe=300.0        # float !
soil_theta_prescribe=1.0      # float !
longitude=0.0                   # float !
latitude=50.0                   # float !

########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere ##
cd "$( dirname "$0" )"


## deal with options and arguments passed to the script ##
while getopts ":p:m:y:o:n:h" opt; do
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
    echo " -p set the PFT (plant_functional_type_id)"
    echo "    Arguments to -p:"
    echo "       any integer value 1-8"
    echo "       note: setPFTAutomatically is set to: no"
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
    echo "       pre: prescribed"
    echo "       std: standard"
    echo "       opt: optimal"
    echo "       res: resistance"
    exit 1
    ;;
  esac
done

## 
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
validateNUptakeScheme="prescribed standard optimal resistance"
validatePUptakeScheme="prescribed standard"
validateCanopyLayerScheme="standard fapar"
validateSoilMoistureScheme="prescribed bucket layer"
validateFlags="F T"
validateOutputInterval="timestep daily weekly monthly yearly"
validatePlantFunctionalTypeID="1 2 3 4 5 6 7 8"
validateOptionForcing="constant generate readfile"
validateSetPFTAutomatically="yes no"
validateSiteDataSet="f c n"

## validate user-defined model-options ##
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"
validateUserConfigValues "$n_uptake_scheme" "$validateNUptakeScheme" "n_uptake_scheme"
validateUserConfigValues "$p_uptake_scheme" "$validatePUptakeScheme" "p_uptake_scheme"
validateUserConfigValues "$canopy_layer_scheme" "$validateCanopyLayerScheme" "canopy_layer_scheme"
validateUserConfigValues "$soil_moisture_scheme" "$validateSoilMoistureScheme" "soil_moisture_scheme"
validateUserConfigValues "$flag_calc_soil_energy" "$validateFlags" "flag_calc_soil_energy"
validateUserConfigValues "$flag_tile_veg_dynamics" "$validateFlags" "flag_tile_veg_dynamics"
validateUserConfigValues "$flag_t_resp_acclimation" "$validateFlags" "flag_t_resp_acclimation"
validateUserConfigValues "$flag_t_jmax_acclimation" "$validateFlags" "flag_t_jmax_acclimation"
validateUserConfigValues "$flag_optimal_Nleaf" "$validateFlags" "flag_optimal_Nleaf"
validateUserConfigValues "$flag_heterogeneous_canopy" "$validateFlags" "flag_heterogeneous_canopy"
validateUserConfigValues "$flag_optimal_allocation" "$validateFlags" "flag_optimal_allocation"
validateUserConfigValues "$flag_optimal_Nfraction" "$validateFlags" "flag_optimal_Nfraction"
validateUserConfigValues "$output_interval_pool" "$validateOutputInterval" "output_interval_pool"
validateUserConfigValues "$output_interval_flux" "$validateOutputInterval" "output_interval_flux"
validateUserConfigValues "$setPFTAutomatically" "$validateSetPFTAutomatically" "setPFTAutomatically"
validateUserConfigValues "$plant_functional_type_id" "$validatePlantFunctionalTypeID" "plant_functional_type_id"
validateUserConfigValues "$flag_oneYearForcing" "$validateFlags" "flag_oneYearForcing"
validateUserConfigValues "$siteDataSet" "$validateSiteDataSet" "siteDataSet"
# not validated: simulation_length_number, output_start_first_day_year, output_end_last_day_year, soil_awc_prescribe, soil_theta_prescribe, latitude, longitude 


## definition of some basic variables
# model forcing - input files for the model
forceSoilLitterFile="soil_forcing_litterfall.dat"
forceSoilPhysicFile="soil_forcing_soilphysics.dat"
forceLandMetFile="climate.dat"

## test if forcing data exist in testbed/ or in the fluxnet/cruncep site-list-folder [if model is not test_canopy] ##
if [ $slmModel != "test_canopy" ] # copy a climate.dat file with PFT "BNE" for all models but test_canopy
  then
    if [ "$siteDataSet" == "n" ]  # if forcing data are neither Fluxnet nor Cruncep dataset
    then  # local forcing file
      if [ $slmModel == "soil" ]
        then
          if [ ! -e ./"$forceSoilLitterFile" ] # if file not exists
            then
              echo "EXIT!  the test_bed/""$forceSoilLitterFile"" does not exist, but the option siteDataSet is set to n !"
              exit
          fi
          if [ ! -e ./"$forceSoilPhysicFile" ] # if file not exists
            then
              echo "EXIT!  the test_bed/""$forceSoilPhysicFile"" does not exist, but the option siteDataSet is set to n !"
              exit
          fi
        else # land plant canopy
          if [ ! -e ./"$forceLandMetFile" ] # if file not exists
            then
              echo "EXIT!  the test_bed/""$forceLandMetFile"" does not exist, but the option siteDataSet is set to n !"
              exit
          fi
      fi 
      setPFTAutomatically="no"
      echo -e '\n' "Info: the option setPFTAutomatically is used with 'no', as you are using the test_bed/climate.dat, which does not provide PFT information. "
    else  # fluxnet/cruncep forcing file
      if [ $siteDataSet == "f" ]
        then
          siteDataSetFull="FLUXNET"
        else
          siteDataSetFull="CRUNCEPv5"
      fi 

      if [ $slmModel == "soil" ]
        then
          forcing="soil_forcing"
        else # land plant canopy
          forcing="met_forcing"
      fi 
      siteFolderBasic="/Net/Groups/BGI/data/Model/SLM1/input/point/"
      siteFolder="$siteFolderBasic""$siteDataSetFull"
      siteFolderForcing="$siteFolder""/""$forcing"
      if ! ls -la $siteFolderForcing | grep $siteName > /dev/null  # if forcing file not found
        then
          echo "EXIT! The site ""$siteName"" you selected for forcing input does not exist in the site data-set folder !"
          exit
      fi
  fi  # ["$siteDataSet" == "n"]
else  # in case: $slmModel == "test_canopy"
  setPFTAutomatically="no"
fi    # [$slmModel != "test_canopy"]
[ $slmModel == "canopy" ] && simulation_length_number=10


## create folder, and copy R scripts, namelist, parameter-list and climate.dat file there ##
mkdir -p ./"$slmModel"
cp ../r_script/postprocessing_"$slmModel".R ./"$slmModel"
cp ../namelist/namelist.slm ./"$slmModel"
if [ $slmModel != "test_canopy" ] # copy a climate.dat file for all models but test_canopy
  then
    if [ "$siteDataSet" == "n" ]
      then
        if [ $slmModel == "soil" ]
          then
            cp ./"$forceSoilLitterFile" ./"$slmModel"/"$forceSoilLitterFile"
            cp ./"$forceSoilPhysicFile" ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp ./"$forceLandMetFile" ./"$slmModel"/"$forceLandMetFile"
        fi 
      else # $siteDataSet == f / c
        if [ $slmModel == "soil" ]
          then
            cp "$siteFolderForcing"/"$siteName"_lit.dat ./"$slmModel"/"$forceSoilLitterFile"
            cp "$siteFolderForcing"/"$siteName"_sph.dat ./"$slmModel"/"$forceSoilPhysicFile"
          else # land plant canopy
            cp "$siteFolderForcing"/"$siteName".dat     ./"$slmModel"/"$forceLandMetFile"
        fi 
    fi
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
sync
return 0                                                                              # get return value with: return_value=$?
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
changeNamelistValue simulation_length_number "$simulation_length_number"
changeNamelistValue output_interval_pool "$output_interval_pool"
changeNamelistValue output_interval_flux "$output_interval_flux"
changeNamelistValue output_start_first_day_year "$output_start_first_day_year"
changeNamelistValue output_end_last_day_year "$output_end_last_day_year"
changeNamelistValue soil_awc_prescribe "$soil_awc_prescribe"
changeNamelistValue soil_theta_prescribe "$soil_theta_prescribe"
changeNamelistValue soil_depth "$soil_depth"
changeNamelistValue set_parameter_values_from_file ".FALSE."   # set to default FALSE because parameter update is not applicable here in test_bed
# get the PFT & lon & lat according to the site selected (reading the all_site_list.dat)
if [ "$setPFTAutomatically" == "yes" ]
  then
    allSiteDatFile="$siteFolder"/all_site_list.dat
    if [ $siteDataSet == "f" ]
      then  # fluxnet sites
        siteNameShort=$(echo "$siteName" | cut -c1-5)             # shortens file name provided, e.g.: CA-Man_1999-2003 -> CA-Man
      else  # cruncep sites
        siteNameShort=$(echo "$siteName" | rev | cut -c11- | rev) # shortens file name provided, e.g.: DOO_00_1984-2013 -> DOO_00
    fi 
    pftName=$( grep "$siteNameShort" "$allSiteDatFile" | awk '{print $4}')
    longitude=$( grep "$siteNameShort" "$allSiteDatFile" | awk '{print $2}')
    latitude=$( grep "$siteNameShort" "$allSiteDatFile" | awk '{print $3}')
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
    changeNamelistValue plant_functional_type_id "$pftID"

  else  # if 'setPFTAutomatically == no' or 'siteDataSet == n'
    changeNamelistValue plant_functional_type_id "$plant_functional_type_id"
fi
changeNamelistValue longitude "$longitude"
changeNamelistValue latitude "$latitude"


## confirmation of finishing this script & starting SLM ##
echo " Nice idea !" #OK"
echo " Let's go with SLM ""$slmModel" #Starting SLM Simulation: ""$slmModel"
echo -n "  ...  "

## run the SLM binary (always compiled as land.x) ##
compiler=$(sed -n 2p ./../../../util/Make.inc.local | awk '{print $3}')
./../../../"$compiler"/bin/land.x 
echo "Simulation finished"
sync

## R script for postprocessing ##
postProcessCanopy=F
echo " Starting postprocessing"
if [ $slmModel == "canopy" ]
 then
  if [ $postProcessCanopy == F ]
   then
    echo "  Oh snap, postprocessing was disabled for SLM CANOPY in the run-slm_test-bed.sh [postProcessCanopy=F]"
  else
    Rscript ./postprocessing_"$slmModel".R
  fi
else
 Rscript ./postprocessing_"$slmModel".R
fi
sync
echo " ... done !"

## end of this script, confirmation about successful simulation is provided by the R script itself ##
