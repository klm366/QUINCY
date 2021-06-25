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
# run this script by executing (at io1): ./compare_fluxnet_diagnostics_single_site.sh.sh
# 
# this script calls the R script ./create_fluxnet_comparison_diagnostics.R 
#  for one particular site of the reference simulation-experiment that is also present in the new simulation experiment
#
# the comparison figure is then available in the "postprocessing/" directory of this site in the new simulation-experiment directory
#
# the new output file fluxnet_diagnostics_comp_[sitename].pdf is created with figures that contain:
#  baseline simulations | new/testing simulations | FLUXNET2 data
#

# get username for the 'BGI/scratch/' directory path to the simulation experiment
username=$(echo $USER)

# provide the location of any previous (e.g. reference/baseline) & your new simulation directories
# note: provide directory paths without a trailing "/" but simply "/path/exp_ID/site-name"
# the reference/baseline simulation-experiment
baselineSiteFolder="/Net/Groups/BGI/work_2/QUINCY/SLM/reference_runs/fluxnet2_siteset/exp_82_baseline_fluxnet_stdBNF_P_rxxxx/DE-Hai"
#baselineSiteFolder="/Net/Groups/BGI/scratch/"$username"/slm/exp_ID/DE-Hai"

# your new simulation experiment to test
newSimulationSiteFolder="/Net/Groups/BGI/scratch/"$username"/slm/exp_ID/DE-Hai"

## modify the location of the all_sites_list.dat if necessary (e.g. when working with a local copy)
## or empty the path to  allSiteList=""  and provide start and last year of fluxnet data for this site manually
## 
## default: allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"
allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"
#allSiteList=""
startYear=2000    # only used if: allSiteList=""
lastYear=2006     # only used if: allSiteList=""


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

# get site name
siteName=$(basename "$newSimulationSiteFolder")

# get start and last year of data from all all_sites_list.dat
if [[ ! -z "$allSiteList" ]]; then   # if $allSiteList is not empty
  startYear=$(grep "$siteName" "$allSiteList" | awk '{print $5}')
  lastYear=$(grep "$siteName" "$allSiteList" | awk '{print $6}')
  sleep 0.4
fi

# start the R script
echo " "
echo ">>  creating comparision plot for ""$siteName"" with fluxnet data from ""$startYear"" to ""$lastYear"
Rscript --vanilla ./create_fluxnet_comparison_diagnostics.R "$startYear" "$lastYear" "$baselineSiteFolder" "$newSimulationSiteFolder" "$siteName"

echo ">>  comparison of Fluxnet diagnostics complete for the site:"
echo "    ""$newSimulationSiteFolder"
