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


#
# run this script by executing (at io1): ./compare_fluxnet_diagnostics_multiple_sites.sh
# 
# this script calls the R script ./create_fluxnet_comparison_diagnostics.R 
#  for all sites of the reference simulation-experiment that are also present in the new simulation experiment
#
# the comparison figures are then available in the directory "x_fluxnet_comparison/" of the new simulation-experiment directory
#
# per site the new output file fluxnet_diagnostics_comp_[sitename].pdf is created with figures that contain:
#  baseline simulations | new/testing simulations | FLUXNET2 data
#

# provide the location of any previous (e.g. reference/baseline) & your new simulation directories
# note: provide directory paths without a trailing "/" but simply "/path/exp_ID"
# the reference simulation-experiment
baselineExpFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/exp_82_baseline_fluxnet_stdBNF_P_rxxxx"
# your new simulation experiment to test
newSimulationExpFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/exp_82_baseline_fluxnet_stdBNF_P_rxxxx"

## modify the location of the all_sites_list.dat if necessary (e.g. when working with a local copy)
## default: allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"
allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"


## set baseline default folder
if [[ "$baselineExpFolder" == "default" ]]; then
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
    if [[ $siteName != "x_binary" ]] && [[ $siteName != "x_fluxnet_comparison" ]] && [[ $siteName != "x_compare"* ]]; then
      
      # if this site is also available in the $newSimulationExpFolder
      if [[ -d "$newSimulationExpFolder""/""$siteName" ]]; then
      
        # get start and last year of data from all all_sites_list.dat
        startYear=$(grep "$siteName" "$allSiteList" | awk '{print $5}')
        lastYear=$(grep "$siteName" "$allSiteList" | awk '{print $6}')
        sleep 0.4

        # start the R script
        echo " "
        echo ">>  creating comparision plot for ""$siteName"" with FLUXNET2 data from ""$startYear"" to ""$lastYear"
        Rscript --vanilla ./create_fluxnet_comparison_diagnostics.R "$startYear" "$lastYear" "$baselineExpFolder""/""$siteName" "$newSimulationExpFolder""/""$siteName" "$siteName"
        echo ">>  comparison of fluxnet diagnostics complete"

        # copy a file of the site comparison fluxnet diagnostics
        cp "$newSimulationExpFolder""/""$siteName""/postprocessing/fluxnet_diagnostics_comp"* "$newSimulationExpFolder""/x_fluxnet_comparison/"
      fi
    fi
done

sync
echo " "
echo " all sites completed"
echo " "
echo " the comparison figures are available in the directory:"
echo "  ""$newSimulationExpFolder""/x_fluxnet_comparison/"
echo " "

## in addition make all files writeable for the group
echo " setting group write-rights to the directories & files within  work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/ "
chmod -R g+w /Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/*
echo " done"
echo " "
