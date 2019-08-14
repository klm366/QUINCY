#!/bin/bash

##    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
##    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
##
##    This file is part of the QUINCY model.
##    The QUINCY model is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    The QUINCY model is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.

# 
# this script runs the R script ./create_fluxnet_comparison_diagnostics.R for simulation results of an individual site
#  essentially creating the new output file FLUXNET_diagnostics_comp_[sitename].pdf with graphics output for:
#  baseline simulations | new/testing simulations | fluxnet data
#


## provide the location of previous/baseline & new/testing simulation results
## note: no trailing "/" but simply "/path/siteName" e.g. "/path/DE-Hai"
baselineSiteFolder="/Net/Groups/BGI/scratch/jengel/slm/exp_r998_fluxdiagn01/DE-Hai"
newSimulationSiteFolder="/Net/Groups/BGI/scratch/jengel/slm/exp_r998_fluxdiagn02/DE-Hai"

## modify the location of the all_sites_list.dat if necessary (e.g. when working with a local copy)
## or empty the path to  allSiteList=""  and provide start and last year of fluxnet data for this site manually
## 
## default: allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"
allSiteList="/Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET2/flux_all_sites_1995-2006_list.dat"
#allSiteList=""
startYear=2000    # only used if allSiteList=""
lastYear=2006     # only used if allSiteList=""



########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

# get site name
siteName=$(basename "$newSimulationSiteFolder")

# get start and last year of data from all all_sites_list.dat
if [ ! -z "$allSiteList" ]   # if $allSiteList is not empty
 then
  startYear=$(grep "$siteName" "$allSiteList" | awk '{print $5}')
  lastYear=$(grep "$siteName" "$allSiteList" | awk '{print $6}')
  sleep 0.4
fi

# start the R script
echo " "
echo ">>  creating comparision plot for ""$siteName"" with fluxnet data from ""$startYear"" to ""$lastYear"
Rscript --vanilla ./create_fluxnet_comparison_diagnostics.R "$startYear" "$lastYear" "$baselineSiteFolder" "$newSimulationSiteFolder" "$siteName"

echo ">>  comparison of Fluxnet diagnostics complete"
