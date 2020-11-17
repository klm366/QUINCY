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

## 
## script to cleanup temporary files and folders after a whole job-array of slm simulations was completed
## 
## this script is copied to /Net/Groups/BGI/scratch/""$username""/slm/$jobName/default_files/ by the prep_run_slm.sh
##  and the bsub runs it from there
## 


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################

## arguments passed to cleanup_slm_hpc_run.sh ##
jobName=$1                        # a user specific alphanumeric string
array=$2                          # was an job array or just one job simulated  T/F
# $3 == singleSiteName provided only if just one single site was run; see below


## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"
slmModel=$(sed -n 3p "$simFold""default_files/config.txt")                                      # QUINCY model: land plant soil canopy
fluxnetTypeTransientTimestepOutput=$(sed -n 6p "$simFold""default_files/config.txt")            # T/F
paramSensitivitySim=$(sed -n 9p "$simFold""default_files/config.txt" | awk '{print $1}')        # T/F  -> parameter sensitivity
compareSimCRUNCEPvsGFDB=$(sed -n 10p "$simFold""default_files/config.txt" | awk '{print $1}')   # T/F  -> compare with global forest data base
createSoilForcFlag=$(sed -n 13p "$simFold""default_files/config.txt")                           # T/F

## move job specific error & out files into the particular job folder ##
if [ "$array" == "T" ]
 # $array == T
 then
  # parameter sensitivity simulations
  if [ $paramSensitivitySim == T ]
    then
      simNumber=$(ls -l "$simFold"error | grep err | wc -l) # number of simulations and number of error files
      # folders start with 1001 for err.2.txt
      for ((i = 0 ; i < $simNumber ; i++))
        do
          errFileCount=$(($i +2))
          simCount=$((1000 + $i + 1))
          mv "$simFold""error/err.""$errFileCount"".txt"  "$simFold""$simCount"_*/error/
          mv "$simFold""out/out.""$errFileCount"".txt"    "$simFold""$simCount"_*/std_out/
      done
  # normal simulations
  else
    # check whether normal or custom site-list was used
    if [ -e "$simFold""default_files/my_site_list.dat" ] # custom site-list exists
     then
      siteList="$simFold""default_files/my_site_list.dat"
     else
      siteList=$(sed -n 1p "$simFold""default_files/config.txt")
    fi
    # move error & out files
    for i in "$simFold"error/*
     do
      siteNumber=$(awk -F "." '{print $2}' <<< $i)
      siteName=$(sed -n "$siteNumber"p "$siteList" | awk '{print $1}')
      mv "$simFold""error/err.""$siteNumber"".txt" "$simFold""$siteName""/error/"
      mv "$simFold""out/out.""$siteNumber"".txt"   "$simFold""$siteName""/std_out/"
    done
  fi
 # $array == no
 else
  singleSiteName=$3               # name of the site
  mv "$simFold""error/err.txt" "$simFold""$singleSiteName""/error/"
  mv "$simFold""out/out.txt"   "$simFold""$singleSiteName""/std_out/"
fi

## create an md5sum file for all output files, i.e. files in the folder "output"
# feature disabled by default, as it takes quite some time for large files, enable manually if needed
# md5sum "$simFold"*/output/* > "$simFold"md5sum-output-files.txt

## save the binary and file with revision number and the parameter table
mkdir -p "$simFold"x_binary
mv "$simFold"default_files/*.x                  "$simFold"x_binary/
mv "$simFold"default_files/revision_number.txt  "$simFold"x_binary/
if [ $paramSensitivitySim == T ]; then
  mkdir -p "$simFold"parameter_list
  mv "$simFold"parameter_sensitivity_files/parameter_table.txt              "$simFold"parameter_list
  mv "$simFold"parameter_sensitivity_files/create_parameter_summary_figs.R  "$simFold"parameter_list
fi

## comparing simulated CRUNCEP data with global forest data base data at daily basis
if [ $compareSimCRUNCEPvsGFDB == T ]; then
  echo " starting sim_cruncep_vs_gfdb-process_daily.R"
  Rscript --vanilla "$simFold"default_files/sim_cruncep_vs_gfdb-process_daily.R $simFold    # creates two .csv files
  sync
  echo " starting sim_cruncep_vs_gfdb-create_figs.R (cruncep CUE sites)"
  Rscript --vanilla "$simFold"default_files/sim_cruncep_vs_gfdb-create_figs.R $simFold      # creates a .pdf with multiple figures
  sync
  echo " starting sim_cruncep_all_vs_gfdb-create_figs.R (all cruncep sites)"
  Rscript --vanilla "$simFold"default_files/sim_cruncep_all_vs_gfdb-create_figs.R $simFold  # creates a .pdf with multiple figures for GPP & NPP
  sync
  echo " finished cruncep vs gfdb data processing and figure creation"
fi

## postprocess timestep output from transient simulations with a fluxnet-type site-set
#  creates: 
#   slm_site_mean.csv & slm_site_sd.csv
#   slm_site_barplots.pdf & slm_site_mean_in_gramm_barplots.csv
if [ "$fluxnetTypeTransientTimestepOutput" == T ]; then
  echo " sim_fluxnet_diagnostics-process_timestep.R"
  Rscript --vanilla "$simFold"default_files/sim_fluxnet_diagnostics-process_timestep.R $simFold  # creates two .csv files
  sleep 2
  echo " finished barplot figure creation of the fluxnet diagnostics"
fi

## create one .pdf file with figs from a particular output file for all sites, for example, vegpoolC
cp "$simFold"default_files/postprocessing_site_overview.R  "$simFold"
cd "$simFold"
Rscript --vanilla ./postprocessing_site_overview.R
sleep 2


## create one .pdf file from the postprocessing .pdf of all sites ##
# currently only for the *plots.pdf but not the *plots_lay.pdf
# do this only if number of folders (sites/simulations) is not more than 100, otherwise 'pdftk' is likely to crash and so this cleanup-script
cd "$simFold"
nSimDir=$(ls -ld */ | wc -l)
if [ $nSimDir -lt 105 ]; then
  mkdir -p ./pdftmp
  count=100
  
  for site in $(ls -d ./*/)
   do
    siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67
    if [ $siteBase != "x_binary" ] && [ $siteBase != "pdftmp" ] && [ $siteBase != "parameter_list" ] && [ $siteBase != "default_files" ] && [ $siteBase != "error" ] && [ $siteBase != "out" ]
    then
      ((count++))
      cp ./"$siteBase"/postprocessing/*plots.pdf  ./pdftmp/"$count"b.pdf
      echo "$jobName"        > ./pdftmp/temp.txt
      echo " "              >> ./pdftmp/temp.txt
      echo " "              >> ./pdftmp/temp.txt
      echo " "              >> ./pdftmp/temp.txt
      echo "~color{1 0 0}"  >> ./pdftmp/temp.txt     # https://stackoverflow.com/questions/18422738/how-to-change-text-color-in-output-from-enscript
      echo $siteBase        >> ./pdftmp/temp.txt
      sleep 1
      enscript -e~ -q -B -f Courier20 -p ./pdftmp/"$count"a.ps ./pdftmp/temp.txt
      sleep 1
      ps2pdf ./pdftmp/"$count"a.ps ./pdftmp/"$count"a.pdf
    fi
  done
  # TODO replace pdftk by qpdf
  pdftk ./pdftmp/*.pdf output ./postprocessing_all_sites.pdf
  rm -r ./pdftmp
fi


## remove some folders not needed anymore ##
sync
rm -r "$simFold"default_files        # includes this cleanup script
rm -r "$simFold"error
rm -r "$simFold"out
if [ $paramSensitivitySim == T ]; then
  rm -r "$simFold"parameter_sensitivity_files
fi


# end cleanup #
