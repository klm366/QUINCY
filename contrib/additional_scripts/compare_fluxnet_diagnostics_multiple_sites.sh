#!/bin/bash

# 
# this script runs the R script ./create_fluxnet_comparison_diagnostics.R for simulation results of multiple sites of one simulation experiment
#  essentially creating per site the new output file FLUXNET_diagnostics_comp_[sitename].pdf with graphics output for:
#  baseline simulations | new/testing simulations | fluxnet data
#


## provide the location of baseline & new/testing simulation folders
## the sites in the new/testing simulation folder can be a subset of the sites in the baseline folder
## you can set baselineExpFolder="default" to use /Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/fluxnet_diagnostics_default
##   which points always to the latest usecase 80 simulations
## 
## note: no trailing "/" but simply "/path/exp_ID" e.g. "/path/exp_fluxnet_test01"
baselineExpFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/FLUXNET/exp_80_baseline_fluxnet_r1130"
newSimulationExpFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/FLUXNET/exp_80_baseline_fluxnet_r1185"

## modify the location of the all_sites_list.dat if necessary (e.g. when working with a local copy)
## default: allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/all_sites_list.dat"
allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1996-2008_list.dat"



########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"


## set baseline default folder
if [ "$baselineExpFolder" == "default" ]
  then
    baselineExpFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/fluxnet_diagnostics_default"
fi


## create a new folder in $newSimulationExpFolder that stores a copy of each new file created
##  and create a info file within this folder
mkdir -p "$newSimulationExpFolder""/x_fluxnet_comparison"
echo "comparison done: "                                     > "$newSimulationExpFolder""/x_fluxnet_comparison/info.txt"
date +%F\ %X                                                >> "$newSimulationExpFolder""/x_fluxnet_comparison/info.txt"
echo " "                                                    >> "$newSimulationExpFolder""/x_fluxnet_comparison/info.txt"
echo "baselineExpFolder: ""$baselineExpFolder"              >> "$newSimulationExpFolder""/x_fluxnet_comparison/info.txt"
echo "newSimulationExpFolder: ""$newSimulationExpFolder"    >> "$newSimulationExpFolder""/x_fluxnet_comparison/info.txt"


## loop over all sites in the $baselineExpFolder 
for site in $(ls -d "$baselineExpFolder"/*/)
  do
    # get the folder name
    siteName=$(basename $site)  # just: DE-Hai or IT-Lav or US-Blo ...

    # if it is a valid site folder
    if [ $siteName != "x_binary" ] && [ $siteName != "x_fluxnet_comparison" ]
      then
        
        # if this site is also available in the $newSimulationExpFolder
        if [ -d "$newSimulationExpFolder""/""$siteName" ]
          then
        
            # get start and last year of data from all all_sites_list.dat
            startYear=$(grep "$siteName" "$allSiteList" | awk '{print $5}')
            lastYear=$(grep "$siteName" "$allSiteList" | awk '{print $6}')
            sleep 0.4

            # start the R script
            echo " "
            echo ">>  creating comparision plot for ""$siteName"" with fluxnet data from ""$startYear"" to ""$lastYear"
            Rscript --vanilla ./create_fluxnet_comparison_diagnostics.R "$startYear" "$lastYear" "$baselineExpFolder""/""$siteName" "$newSimulationExpFolder""/""$siteName" "$siteName"
            echo ">>  comparison of Fluxnet diagnostics complete"

            # copy a file of the site comparison fluxnet diagnostics
            cp "$newSimulationExpFolder""/""$siteName""/postprocessing/fluxnet_diagnostics_comp"* "$newSimulationExpFolder""/x_fluxnet_comparison/"
        fi
    fi
done

sync
echo " "
echo " all sites completed"
echo " "

chmod g+w -R /Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/FLUXNET/*
echo " and also set group write-rights to the folders within  work_2/QUINCY/SLM/ReferenceRuns/FLUXNET/ "
echo " "
