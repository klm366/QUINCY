#!/bin/bash

##
## create comparison figures from default output of two simulation experiments
##
## HowTo:
##  copy this bash script and the 'compare_two_quincy_exp_create_figs.r' into a new directory at your scratch (i.e. BGI/scratch/USER/)
##  adapt the below four settings to your needs
##  execute this bash script at either io1 or macadamia
##  the comparison figures may be created in the 'figs' directory (the 'old-output' and 'new-output' are just temporary directories)
##

## get username (for default quincy output directory at the BGI scratch)
username=$(echo $USER)

###################################################################################################
## modify settings 1 to 4 as you see fit
###################################################################################################

# 1 #  output-period of output files used for comparison figures
#outF="yearly"
#outF="monthly"
#outF="weekly"
#outF="daily"
outF="timestep"

# 2 #  set IDs for old and new experiment (used in R script for the figure legends)
# "old" is the referenence you want to compare your "new" Code with
# usually these IDs would contain the revision and a branch-related word (no spaces!)
# e.g.: oldID="r2xxx_trunk"
oldID=""
newID=""

# 3 #  experiment-directory names
# e.g.: old="exp_82_baseline_fluxnet_stdBNF_P_r2xxx"
old=""
new=""

# 4 #  set path to experiments
# old simulation
#baseDirO="/Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/"
baseDirO="/Net/Groups/BGI/scratch/""$username""/slm/"

# new simulation
#baseDirN="/Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/"
baseDirN="/Net/Groups/BGI/scratch/""$username""/slm/"

###################################################################################################
## do not change any below this line unless you know what you are doing
###################################################################################################


###################################################################################################
## some automatic settings and tests, and feedback about settings
###################################################################################################

## change into the directory where this script is located, and set as working directory
cd "$( dirname "$0" )"
workDir=$(pwd)
cd "$workDir"

## this script may run only at io1 or macadamia
node=$(hostname -s)
if [ "$node" != "io1" ] && [ "$node" != "macadamia" ]; then
  echo " EXIT! run this script at either io1 or macadamia"
  echo " You are at node: ""$node"
  exit 1
fi

## R script for figure creation - define filename and test presence
createFigsScriptR="compare_two_quincy_exp_create_figs.r"
if [ ! -e "$createFigsScriptR" ]; then   # if file is NOT present
  echo " EXIT! the R script for creating the figures is not present in this directory"
  echo " copy the script '""$createFigsScriptR""' from the repo and go again"
  exit 1
fi

## feedback
echo "working directory: ""$workDir"
echo "using output files with output-period: ""$outF"
echo "old experiment with ID ""$oldID"":"
echo "  ""$baseDirO""$old"
echo "new experiment with ID ""$newID"":"
echo "  ""$baseDirN""$new"
echo " "


###################################################################################################
## remove old output and figure files (& create directories if not present already)
###################################################################################################
mkdir -p ./old-output
mkdir -p ./new-output
mkdir -p ./figs
rm -f ./old-output/*
rm -f ./new-output/*
rm -f ./figs/*
sync


###################################################################################################
## loop over sites of the new simulation that have a matching site in the old simulation
###################################################################################################
for site in $(ls -d "$baseDirN""$new"/*/); do
  ## get folder/site name
  siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67
  ## exclude non-site folders
  ## and run only if there is a matching site (directory) in the "old simulation" experiment directory
  if [ $siteBase != "x_binary" ] && [ $siteBase != "x_fluxnet_comparison" ] && [ $siteBase != "0_compare" ] && [ -d "$baseDirO""$old"/"$siteBase" ]; then
    # test if vegfluxC_[output-period].txt is present in the new simulation
    if [ ! -e "$baseDirN""$new""/""$siteBase""/output/"vegfluxC*"$outF"".txt" ]; then   # if file is NOT present
      echo " EXIT! The following file is not present:"
      echo " ""$baseDirN""$new""/""$siteBase""/output/"vegfluxC*"$outF"".txt"
      echo " maybe because the output-period you selected does not match the simulations output."
      exit 1
    fi
    # copy all output files but not the "soil_layer" files (and neither the veg_canopy that has a custom internal structure)
    rsync --exclude="*_lay_*" --exclude="*veg_canopy*" "$baseDirO""$old""/""$siteBase"/output/*"$outF".txt    ./old-output/
    rsync --exclude="*_lay_*" --exclude="*veg_canopy*" "$baseDirN""$new""/""$siteBase"/output/*"$outF".txt    ./new-output/
    sync

    # run R script to create the figures
    # if old/new ID are empty, set default
    if [ -z "$oldID" ]; then
      oldID="old"
    fi
    if [ -z "$newID" ]; then
      newID="new"
    fi
    Rscript --vanilla "$createFigsScriptR" "$workDir" "$oldID" "$newID" "$outF"
    sleep 0.2
    sync
    rm -f ./old-output/*
    rm -f ./new-output/*

    # move created .pdf files into the ./figs directory
    compFlag="F"
    diffFlag="F" 
    if [ -e ./figs/comp_results_two_exp.pdf ]; then
      mv ./figs/comp_results_two_exp.pdf  ./figs/"$siteBase"_compare.pdf
      compFlag="T"
    fi
    if [ -e ./figs/diff-prop_two-exp_new-minus-old.pdf ]; then
      mv ./figs/diff-prop_two-exp_new-minus-old.pdf  ./figs/"$siteBase"_diff-prop_new-minus-old.pdf
      diffFlag="T"
    fi
    echo "done: ""$siteBase"
  fi
done

## finish
rm -r ./old-output/
rm -r ./new-output/
sync
echo " "
echo " finished !"

