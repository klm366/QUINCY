#!/bin/bash

## 20160810
## Where: you can start this script from anywhere (no need to cd into a specific directory), when logged into a cluster node such as pc028
## What:  this script starts one of the run_xx.R scripts to run a specific SLM model, e.g. land plant canopy
## HowTo: either specify the model a) here within this bash script or b) pass is as argument to this script
## a) script option slmModel
## b) run: run-slm_test-bed.sh modelName
##   b has priority over a


########################################################################################################
####### simulation specific options you may want to change #############################################
########################################################################################################

slmModel="plant"   # specific SLM model [default: "land"]; possible values: land plant canopy test_canopy




########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere ##
cd "$( dirname "$0" )"


## if slmModel is passed as argument to this script, use this over the above defined value ##
if [ ! -z "$1" ] # if $1 is not empty -> argument passed to this script
  then
    slmModel=$1
fi


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
validateSlmModel="land plant canopy test_canopy"  # soil 


## validate user-defined model-options ##
validateUserConfigValues "$slmModel" "$validateSlmModel" "slmModel"


## manipulate namelist.slm accordingly to above options ##
 # function: pass the variable_name ($1) & the new_value ($2)
 # this function assumes a specific relative path to the namlist.slm to be valid
function changeNamelistValue {
oldValue=$(grep -w $1 ../namelist/namelist.slm | awk -F "=" '{print $2}' | tr -d '"') # get the old value that should be replaced
lineNum=$(grep -wn $1 ../namelist/namelist.slm | cut -d : -f 1)                       # get line number of variable $1
sed -i "$lineNum"s/$oldValue/$2/ ../namelist/namelist.slm                             # replace the old value in line $lineNum with $2
return 0                                                                              # get return value with: return_value=$?
}
 # examples:
 #changeNamelistValue simulation_length_number 300
 #changeNamelistValue output_interval_flux monthly    # quotation marks not needed
 #changeNamelistValue flag_calc_soil_moisture .FALSE.
 #changeNamelistValue soil_theta_prescribe 122.7

changeNamelistValue slm_model_name "$slmModel"
sync

## confirmation of starting the R script ##
echo " starting the R script now running SLM: ""$slmModel"


## cd into the specific folder and run the R script ##
cd ./"$slmModel"
Rscript ./run_"$slmModel".R



## end of this script, confirmation about successful simulation is provided by the R script itself ##
