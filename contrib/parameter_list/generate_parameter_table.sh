#!/bin/bash

##
## creates a table (param_table.txt) with sets of parameter value combinations used as input for SLM parameter sensitivity simulations 
## see prep_run_slm.sh option: paramSensitivitySim
## 
## input:   slm_parameter_sensitivity_input.txt (content must comply with the Main Task: https://projects.bgc-jena.mpg.de/QUINCY/wiki/SLMParameterSensitivitySimulations)
## output:  param.listD
##          param_table.txt
##
## 1st version: Jan Engel MPI-BGC Jena 20161006
##
##

########################################################################################################
####### options you might want to change ###############################################################
########################################################################################################

# Main Tasks
#   1: you provide parameters and specific values per parameter (up to 4 parameters, number of values can differ between parameters)
#   2: you provide parameters, min and max values and number of values per parameter (same conditions as for 1)
#   3: you provide parameters, min and max values and overall number of samples to be taken, the rest is done by LHS using R
# Note: you can pass the mainTask number as argument to this script; for Main Task 2: ./path-to-script/generate_parameter_table.sh 2
mainTask=1


# Input File
#   defaults to: ./slm_parameter_sensitivity_input.txt
#   leave paramInputFile blank to use the default
#   you can overwrite this default by providing a valid path/file.txt relative to this location contrib/parameter_list/
paramInputFile="" 


# Output Files (no options here)
#   the two files will be saved to this folder (contrib/parameter_list/ ) and deleted during the next run of this script
#   the prep_run_slm.sh expects these files to be located here in case you set the option: paramSensitivitySim="T"
#   filenames:
#     param.listD
#     param_table.txt

########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"
workingDir=$(pwd)       # needed for the below R script

## test for and read argument to this script ##
# this overwrites the above defined value of $mainTask
if [ ! -z "$1" ] # if $1 is not empty -> argument passed to this script
  then
    mainTask=$1
fi

## test and set input file ##
if [ -z "$paramInputFile" ] # if $paramInputFile is empty
  then
    paramInputFile="./slm_parameter_sensitivity_input.txt"
  else
    if [ ! -e "$paramInputFile" ]     # exit if non-existing
      then
    	  echo "EXIT!  I cannot find the custom input file you specified. Please check: ""$paramInputFile"
    	  exit
    fi
fi

## validate main task ##
validate="$(grep -w "$mainTask" <<< "1 2 3")"
if [ -z "${validate// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  echo "EXIT! > ""$mainTask"" <  is not a valid value for parameter mainTask, please change and restart generate_parameter_table.sh"
  exit
 else
  if [ -z $mainTask ]
   then
    echo "EXIT!  The value for parameter mainTask is empty but must contain a valid string, please change and restart generate_parameter_table.sh"
    exit
  fi
fi

## validate param.list ##
if [ ! -e ./param.list ]     # exit if non-existing
  then
    echo "EXIT!  I cannot find the file contrib/parameter_list/param.list."
    exit
fi

## set filenames and variables##
paramList="./param.list"
paramListDefault="./param.listD"
paramTable="./param_table.txt"

paramNumber=0
paramValSets=0

## remove possible white spaces at line start and end, and remove possible empty lines ##
sed -i 's/^ *//'  "$paramInputFile"
sed -i 's/ *$//'  "$paramInputFile"
sed -i '/^\s*$/d' "$paramInputFile"
sync

## validate the parameter in the input file ##
paramListLineNumbers=$(wc -l "$paramInputFile" | awk '{print $1}') # number of parameters in the input file ($paramInputFile)
for (( a=1; a<=$paramListLineNumbers; a++ )) # each line number
  do
    parameterName=$(sed -n "$a"p "$paramInputFile" | awk '{print $1}')  # -n  == quiet/silent
    if ! cat $paramList | grep $parameterName > /dev/null  # if forcing file not found
      then
        echo "EXIT!  The parameter >""$parameterName""< from your ""$paramInputFile"" does not exist in ""$paramList"" !"
        echo " Please add it to both files: ""$paramList"" and mo_shared_set_parameters.f90."
        echo " Ask Jan in case you have further questions or need help with this."
        exit
    fi
done

## remove previous versions of output files ##
rm -f $paramListDefault
rm -f $paramTable
echo "Previous output files have been removed (param.listD & param_table.txt, if existent)."
echo " "

## get the number of lines that is the number of parameters ##
paramNumber=$(wc -l "$paramInputFile" | awk '{print $1}')
declare -a paramValueNumbersA=($paramNumber)

## test and notification of the user about the number of parameters ##
if [ $paramNumber -gt 4 ]
  then
    if [ $mainTask -lt 3 ] # for Main Tasks 1 & 2
      then
        echo "EXIT !  Sorry, the number of parameters is currently limited to 4 for Main Tasks 1 & 2 (you have ""$paramNumber"")."
        exit
    fi
  else
    echo "The number of parameters you want to vary is: ""$paramNumber"
fi

## create an array with parameter names used as header in $paramTable ##
paramTableHead[0]=0
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  paramTableHead[$line]=$(sed -n "$line"p "$paramInputFile" | awk '{print $1}')
done

## function definition ###############################################
## a function that combines some steps of work for Main Tasks 1 & 2 ##
function mainTask1and2 {
## test and notification of the user about the number of variations ##
if [ $paramValSets -gt 994 ]
  then
    echo "EXIT !  You want to do ""$paramValSets"" simulations. That is too much. A number lower than 995 is necessary."
    exit
  else
    echo "The number of combinations of parameter values (i.e. simulations) is: ""$paramValSets"
fi

## create the param.listD file ##
valDefault="NAN"
paramActive="ACTIVE"
cp $paramList $paramListDefault
sync

# mark all used parameters with ACTIVE 
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  paramName=$(sed -n "$line"p "$paramInputFile" | awk '{print $1}')
  sed -i /$paramName/s/$valDefault/"$paramActive"/ $paramListDefault
done

# delete the parameters not used
sed -i "/$valDefault/d" "$paramListDefault"

return 0 
}  # end function ####################################################


###################################################################################################
## Main Tasks 1 ##                                                                               ##
if [ $mainTask == 1 ]; then                                                                      ##
###################################################################################################

## get numbers of values per parameter and cals number of different parameter value sets ##
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  valNum=$(sed -n "$line"p "$paramInputFile" | grep -o " " | wc -l) # count the number of separators between columns
  if [ $i -lt 1 ]
    then
      paramValSets=$valNum
    else
      paramValSets=$(($paramValSets * $valNum))
  fi
done

## call function
mainTask1and2

## create the $paramTable ##
# bash function to combine elements of one string separated by a certain string, such as , or /
function join_by { local IFS="$1"; shift; echo "$*"; } 

if [ $paramNumber == 1 ]
  then
    temp=$(sed -n 1p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param1Val=( $temp )
    joinedParam1Val=$(join_by , ${param1Val[@]})    
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")
paramT<-expand.grid( ${paramTableHead[1]} = c( $joinedParam1Val ))
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 2 ]
  then
    temp=$(sed -n 1p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param1Val=( $temp )
    joinedParam1Val=$(join_by , ${param1Val[@]})    
    temp=$(sed -n 2p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param2Val=( $temp )
    joinedParam2Val=$(join_by , ${param2Val[@]})  
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")
paramT<-expand.grid( ${paramTableHead[1]} = c( $joinedParam1Val ), ${paramTableHead[2]} = c( $joinedParam2Val ))
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 3 ]
  then
    temp=$(sed -n 1p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param1Val=( $temp )
    joinedParam1Val=$(join_by , ${param1Val[@]})    
    temp=$(sed -n 2p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param2Val=( $temp )
    joinedParam2Val=$(join_by , ${param2Val[@]})  
    temp=$(sed -n 3p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param3Val=( $temp )
    joinedParam3Val=$(join_by , ${param3Val[@]})  
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")
paramT<-expand.grid( ${paramTableHead[1]} = c( $joinedParam1Val ), ${paramTableHead[2]} = c( $joinedParam2Val ), ${paramTableHead[3]} = c( $joinedParam3Val ))
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 4 ]
  then
    temp=$(sed -n 1p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param1Val=( $temp )
    joinedParam1Val=$(join_by , ${param1Val[@]})    
    temp=$(sed -n 2p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param2Val=( $temp )
    joinedParam2Val=$(join_by , ${param2Val[@]})  
    temp=$(sed -n 3p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param3Val=( $temp )
    joinedParam3Val=$(join_by , ${param3Val[@]})  
    temp=$(sed -n 4p "$paramInputFile" | awk '{print substr($0, index($0, $2))}')  # do not consider the 1st column
    param4Val=( $temp )
    joinedParam4Val=$(join_by , ${param4Val[@]})  
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")
paramT<-expand.grid( ${paramTableHead[1]} = c( $joinedParam1Val ), ${paramTableHead[2]} = c( $joinedParam2Val ), ${paramTableHead[3]} = c( $joinedParam3Val ), ${paramTableHead[4]} = c( $joinedParam4Val ))
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

# add a zero to the start of the header line 
sed -i '1s/^/0 /' $paramTable


#############################
fi                          #
## End of Main Tasks 1 ######
#############################


###################################################################################################
## Main Tasks 2 ##                                                                               ##
if [ $mainTask == 2 ]; then                                                                      ##
###################################################################################################

## get numbers of values per parameter and cals number of different parameter value sets ##
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  valNum=$(sed -n "$line"p "$paramInputFile" | awk '{print $4}') # get the number of values per parameter from column 4
  if [ $i -lt 1 ]
    then
      paramValSets=$valNum
    else
      paramValSets=$(($paramValSets * $valNum))
  fi
done

## call function
mainTask1and2

## validate basic style of $paramInputFile ##
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  valNum=$(sed -n "$line"p "$paramInputFile" | grep -o " " | wc -l) # count the number of separators between columns
  if [ $valNum != 3 ]
    then
      echo "Something seems wrong with your input file in line number ""$line""."
      echo " You selected Main Task 2, so each line may contain: Variable_name Min_value Max_value Number_of_values (space separated)."
      echo " EXIT!  Please check your input file."
      exit
  fi
done

## create $paramTable
if [ $paramNumber == 1 ]
  then
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")

inputFile<-read.table("$paramInputFile",sep=" ",header=F)
param1<-seq(inputFile[1,2], inputFile[1,3], length.out=inputFile[1,4])

paramT<-expand.grid( ${paramTableHead[1]} = param1)
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 2 ]
  then
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")

inputFile<-read.table("$paramInputFile",sep=" ",header=F)
param1<-seq(inputFile[1,2], inputFile[1,3], length.out=inputFile[1,4])
param2<-seq(inputFile[2,2], inputFile[2,3], length.out=inputFile[2,4])

paramT<-expand.grid( ${paramTableHead[1]} = param1, ${paramTableHead[2]} = param2)
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 3 ]
  then
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")

inputFile<-read.table("$paramInputFile",sep=" ",header=F)
param1<-seq(inputFile[1,2], inputFile[1,3], length.out=inputFile[1,4])
param2<-seq(inputFile[2,2], inputFile[2,3], length.out=inputFile[2,4])
param3<-seq(inputFile[3,2], inputFile[3,3], length.out=inputFile[3,4])

paramT<-expand.grid( ${paramTableHead[1]} = param1, ${paramTableHead[2]} = param2, ${paramTableHead[3]} = param3)
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

if [ $paramNumber == 4 ]
  then
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")

inputFile<-read.table("$paramInputFile",sep=" ",header=F)
param1<-seq(inputFile[1,2], inputFile[1,3], length.out=inputFile[1,4])
param2<-seq(inputFile[2,2], inputFile[2,3], length.out=inputFile[2,4])
param3<-seq(inputFile[3,2], inputFile[3,3], length.out=inputFile[3,4])
param4<-seq(inputFile[4,2], inputFile[4,3], length.out=inputFile[4,4])

paramT<-expand.grid( ${paramTableHead[1]} = param1, ${paramTableHead[2]} = param2, ${paramTableHead[3]} = param3, ${paramTableHead[4]} = param4)
write.table(format(paramT,scientific=F), file="$paramTable", append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)
EOF
# End R script
fi

# add a zero to the start of the header line 
sed -i '1s/^/0 /' $paramTable


#############################
fi                          #
## End of Main Tasks 2 ######
#############################


###################################################################################################
## Main Tasks 3 ##                                                                               ##
if [ $mainTask == 3 ]; then                                                                      ##
###################################################################################################

## 2do !! line 100 -> check for variables validity

echo "Sorry, the latin hypercube sampling (Main Tasks 3) has not been implemented yet."
echo "Please wait for a future release of SLM implementing this feature, or send a reminder to Jan."
echo "EXIT !"
exit


#############################
fi                          #
## End of Main Tasks 3 ######
#############################


echo " "
echo "OK. The script has finished."
echo "  Two new files have been created (param.listD & param_table.txt)."
echo "  You could run the prep_run_slm.sh with paramSensitivitySim == T now."

# end generate_parameter_table.sh # 
