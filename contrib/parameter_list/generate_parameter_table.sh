#!/bin/bash

##
## Info: 
##      * run this script to create the two files needed for starting the parameter sensitivity analysis with the prep_run_slm.sh [option: paramSensitivitySim]
##      * the prep_run_slm.sh expects these files to be located here (contrib/parameter_list/ ) in case you set the option: paramSensitivitySim="T"
##      * this folder (contrib/parameter_list/ ) already contains the two files: generate_parameter_table.sh & param.list
##      * the variables you want to use for your parameter sensitivity analysis must already be present in the param.list, if not, ask Jan or Sönke
##      * you need to create an additional input file (see below for details)
##      * set the mainTask to 1, 2 or 3 (see below for details)
## 
## Input you have to create:   
##      slm_parameter_sensitivity_input.txt 
##        * content of this file must comply with the Main Task: https://projects.bgc-jena.mpg.de/QUINCY/wiki/SLMParameterSensitivitySimulations
##        * you could use a different file name and set the option paramInputFile accordingly
## 
## Output created by this script:  
##      param.listD           
##        * this is a namelist containing the variables you use with value "ACTIVE", is a subset of the variables in param.list
##      param_table.txt
##        * this file contains sets of parameter value combinations used as input for SLM parameter sensitivity simulations
##        * this file is also meant as documentation of the parameter values, a copy is saved in the simulation results folder
##
## 1st version: Jan Engel MPI-BGC Jena 20161006
## 2nd version: Jan Engel MPI-BGC Jena 20180403 add LHS
##

########################################################################################################
####### options you might want to change ###############################################################
########################################################################################################

# Main Tasks
#   1: you provide parameters and specific values per parameter (up to 4 parameters, number of values can differ between parameters)
#   2: you provide parameters, min and max values and number of values per parameter (same conditions as for 1)
#   3: you provide parameters, min and max values and number of parameter sets / simulations, the rest is done by LHS using R
# Note: you can pass the mainTask number as argument to this script; for Main Task 2: ./contrib/generate_parameter_table.sh 2
mainTask=1


# Input file for this script
#   defaults to: ./slm_parameter_sensitivity_input.txt
#   leave paramInputFile blank to use the default
#   you can overwrite this default by providing a valid filename located in this folder, such as, my_param_sens_01.txt
paramInputFile="" 


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"
workingDir=$(pwd)       # needed for the below R script

## test for presence and read argument to this script ##
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
    	  exit 1
    fi
fi

## validate main task ##
validate="$(grep -w "$mainTask" <<< "1 2 3")"
if [ -z "${validate// }" ]  # the search term is available yes / no  [-z test is explicitly for: the length of string is zero]
 then 
  echo "EXIT! > ""$mainTask"" <  is not a valid value for parameter mainTask, please change and restart generate_parameter_table.sh"
  exit 1
else
 if [ -z $mainTask ]        # test if any value is provided to mainTask
  then
   echo "EXIT!  The value for parameter mainTask is empty but must contain a valid string, please change and restart generate_parameter_table.sh"
   exit 1
 fi
fi

## validate existence of param.list ##
if [ ! -e ./param.list ]     # exit if non-existing
  then
    echo "EXIT!  I cannot find the file contrib/parameter_list/param.list."
    exit 1
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

## remove previous versions of output files ##
rm -f $paramListDefault
rm -f $paramTable
echo " "
echo "Previous output files have been removed (param.listD & param_table.txt, if existent)."
echo " "

## get the number of lines that is the number of parameters and create an array with size of number of parameters ##
# in case of mainTask=3 the last line may not contain a parameter but the number of values per parameter
paramNumber=$(wc -l "$paramInputFile" | awk '{print $1}')
if [ "$mainTask" == 3 ]; then paramNumber=$(($paramNumber - 1)); fi
declare -a paramValueNumbersA=($paramNumber)

## validate the parameter names in the input file ##
for (( a=1; a<=$paramNumber; a++ )) # each line number
  do
    parameterName=$(sed -n "$a"p "$paramInputFile" | awk '{print $1}')  # -n  == quiet/silent
    if ! cat $paramList | grep $parameterName > /dev/null  # if parameter name not found in file param.list
      then
        echo "EXIT!  The parameter >""$parameterName""< from your ""$paramInputFile"" does not exist in ""$paramList"" !"
        echo " Please add it to both files: ""$paramList"" and mo_shared_set_parameters.f90."
        echo " Ask Jan or Sönke in case you have further questions or need help with this."
        exit 1
    fi
done

## test and notification of the user about the number of parameters ##
if [ $paramNumber -gt 4 ]
  then
    if [ $mainTask -lt 3 ] # for Main Tasks 1 & 2
      then
        echo "EXIT !  Sorry, the number of parameters is limited to 4 for Main Tasks 1 & 2 (you have ""$paramNumber"")."
        echo "        Please consider using Main Task 3 that applies the LHS and has no such limitations."
        exit 1
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

######################################################################
## function definition ###############################################
## a function that combines some steps of work for Main Tasks 1 2 3 ##
function f_showNumberSetsPrepParamlistD {
## test and notification of the user about the number of variations ##
if [ $paramValSets -gt 994 ]
  then
    echo "EXIT !  You want to run simulations for ""$paramValSets"" parameter sets."
    echo "        That is too much. A number lower than 995 is necessary."
    exit 1
  else
    echo "The number of parameter sets (i.e. simulations) is: ""$paramValSets"
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
######################################################################


###################################################################################################
## Main Task 1 ##                                                                                ##
if [ $mainTask == 1 ]; then                                                                      ##
###################################################################################################

## get numbers of samples per parameter and calc number of parameter value sets ##
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
f_showNumberSetsPrepParamlistD

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
## End of Main Task 1 #######
#############################


###################################################################################################
## Main Tasks 2 ##                                                                               ##
if [ $mainTask == 2 ]; then                                                                      ##
###################################################################################################

## validate basic style of $paramInputFile ##
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  valNum=$(sed -n "$line"p "$paramInputFile" | grep -o " " | wc -l) # count the number of separators between columns
  if [ $valNum != 3 ]
    then
      echo "Something seems wrong with your input file in line number ""$line""."
      echo " You selected Main Task 2, so each line may contain: Variable_name Min_value Max_value Number_of_values (space separated)."
      echo " EXIT!  Please check your input file (""$paramInputFile"")."
      exit 1
  fi
done

## get numbers of samples per parameter and calc number of parameter value sets ##
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
f_showNumberSetsPrepParamlistD


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
## End of Main Task 2 #######
#############################


###################################################################################################
## Main Task 3 ##                                                                                ##
if [ $mainTask == 3 ]; then                                                                      ##
###################################################################################################

## validate basic style of $paramInputFile ##
for ((i = 0 ; i < $paramNumber ; i++))	
do
  line=$(($i + 1))
  valNum=$(sed -n "$line"p "$paramInputFile" | grep -o " " | wc -l) # count the number of separators between columns
  if [ $valNum != 2 ]
    then
      echo "Something seems wrong with your input file in line number ""$line""."
      echo " You selected Main Task 3, so each line may contain: Variable_name Min_value Max_value (space separated)."
      echo " With the last line containing in addition the keyword param_sets and the number of parameter sets / simulations."
      echo " EXIT!  Please check your input file (""$paramInputFile"")."
      exit 1
  fi
done

## validate and get number of samples ##
paramSetsLineNum=$(($paramNumber + 1))
validParamSetsKeyword=$(sed -n "$paramSetsLineNum"p "$paramInputFile" | awk '{print $1}')
if [ "$validParamSetsKeyword" == "param_sets" ]; then
  paramSets=$(sed -n "$paramSetsLineNum"p "$paramInputFile" | awk '{print $2}')
else
  echo "The keyword < param_sets > is missing from your input file."
  echo "You are using Main Tasks 3, where the last line of your input file may contain: param_sets X (X being a number)."
  echo "EXIT !"
  exit 1
fi

## pass the number of parameter sets to the variable used in the function f_showNumberSetsPrepParamlistD ##
paramValSets="$paramSets"

## call function
f_showNumberSetsPrepParamlistD


## create $paramTable
echo " "
echo " Now the values are sampled for each of the parameters using R."
echo " That might take some minutes! depending on the number of parameters and parameter sets."
echo " (For three parameters and 80 parameter sets it takes two and a half minutes at pc027.)"
# R script  # for debug use: R --slave --no-save << EOF
R --slave --silent --no-save << EOF
setwd("$workingDir")
library(lhs)

# get input file
inputFile <- read.csv("$paramInputFile",
                      sep="", header=F, stringsAsFactors=FALSE, as.is=T)
inputFile <- inputFile[1:(nrow(inputFile)-1),]  # remove the last row containing the info on parameter sets

# run LHS
sampleMatrix <- geneticLHS($paramSets,
                           $paramNumber,
                           pop=50, gen=20000, pMut=0.2, criterium="S", verbose=F)

# transform [0,1] values into [min,max] values as defined in the input file
for (p in 1:$paramNumber){
  minVal  <- inputFile[p,2]
  maxVal  <- inputFile[p,3]
  diffVal <- maxVal - minVal
  for (s in 1:$paramSets){
    sampleMatrix[s,p] <- sampleMatrix[s,p] * diffVal + minVal
  }
}

# change matrix to data.frame and write file
sampleMatrix           <- as.data.frame(sampleMatrix)
colnames(sampleMatrix) <- inputFile[,1]

write.table(format(sampleMatrix,scientific=F), file="$paramTable", 
            append=FALSE, sep = " ", eol="\n",quote = FALSE, row.names=TRUE, col.names=TRUE)

EOF
# End R script

# add a zero to the start of the header line 
sed -i '1s/^/0 /' "$paramTable"
sync

#############################
fi                          #
## End of Main Task 3 #######
#############################


echo " "
echo "OK. The script has finished."
echo "  Two new files have been created (param.listD & param_table.txt)."
echo "  You could now run the prep_run_slm.sh with paramSensitivitySim == T."
echo " "

## end generate_parameter_table.sh ##
