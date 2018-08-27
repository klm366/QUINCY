#!/bin/bash

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
slmModel=$(sed -n 3p "$simFold""default_files/config.txt")                                      # SLM model: land plant soil canopy
specialFluxnetOutput=$(sed -n 6p "$simFold""default_files/config.txt")                          # T/F  -> fluxnet diagnostics from timestep output
paramSensitivitySim=$(sed -n 9p "$simFold""default_files/config.txt" | awk '{print $1}')        # T/F  -> parameter sensitivity
compareSimCRUNCEPvsGFDB=$(sed -n 10p "$simFold""default_files/config.txt" | awk '{print $1}')   # T/F  -> compare with global forest data base


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
md5sum "$simFold"*/output/* > "$simFold"md5sum-output-files.txt

## save the binary and the parameter table
mkdir -p "$simFold"x_binary
mv "$simFold"default_files/*.x "$simFold"x_binary
if [ $paramSensitivitySim == T ]; then
  mkdir -p "$simFold"parameter_list
  mv "$simFold"parameter_sensitivity_files/parameter_table.txt              "$simFold"parameter_list
  mv "$simFold"parameter_sensitivity_files/create_parameter_summary_figs.R  "$simFold"parameter_list
fi

## run two R Scripts for comparing simulated CRUNCEP data with global forest data base data at daily basis
if [ $compareSimCRUNCEPvsGFDB == T ]; then
  echo " starting sim_cruncep_vs_gfdb-process_daily.R"
  Rscript --vanilla "$simFold"default_files/sim_cruncep_vs_gfdb-process_daily.R $simFold  # creates two .csv files
  sync
  echo " starting sim_cruncep_vs_gfdb-create_figs.R"
  Rscript --vanilla "$simFold"default_files/sim_cruncep_vs_gfdb-create_figs.R $simFold    # creates a .pdf with multiple figures
  sync
  echo " finished cruncep vs gfdb data processing and figure creation"
fi

## run one R Script for postprocessing the fluxnet diagnostics with barplots
if [ "$specialFluxnetOutput" == T ]; then
  echo " sim_fluxnet_diagnostics-process_timestep.R"
  Rscript --vanilla "$simFold"default_files/sim_fluxnet_diagnostics-process_timestep.R $simFold  # creates two .csv files
  sync
  echo " finished barplot figure creation of the fluxnet diagnostics"
fi


## create one .pdf file with figs from one output file for all sites, for example, vegpoolC
cp "$simFold"default_files/postprocessing_site_overview.R  "$simFold"
cd "$simFold"
Rscript --vanilla ./postprocessing_site_overview.R


## create one .pdf file from the postprocessing .pdf of all sites ##
# currently only for the *plots.pdf but not the *plots_lay.pdf
cd "$simFold"
mkdir -p ./pdftmp
count=100

for site in $(ls -d ./*/)
 do
  siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67
  if [ $siteBase != "x_binary" ] && [ $siteBase != "pdftmp" ] && [ $siteBase != "parameter_list" ]
  then
    ((count++))
    cp ./"$siteBase"/postprocessing/*plots.pdf  ./pdftmp/"$count"b.pdf
    echo "$jobName"        > ./pdftmp/temp.txt
    echo " "              >> ./pdftmp/temp.txt
    echo " "              >> ./pdftmp/temp.txt
    echo " "              >> ./pdftmp/temp.txt
    echo "~color{1 0 0}"  >> ./pdftmp/temp.txt     # https://stackoverflow.com/questions/18422738/how-to-change-text-color-in-output-from-enscript
    echo $siteBase        >> ./pdftmp/temp.txt
    sync
    enscript -e~ -q -B -f Courier20 -p ./pdftmp/"$count"a.ps ./pdftmp/temp.txt
    sync
    ps2pdf ./pdftmp/"$count"a.ps ./pdftmp/"$count"a.pdf
  fi
done
pdftk ./pdftmp/*.pdf output ./postprocessing_all_sites.pdf
rm -r ./pdftmp


## remove some folders not needed anymore ##
rm -r "$simFold"default_files        # includes this cleanup script
rm -r "$simFold"error
rm -r "$simFold"out
if [ $paramSensitivitySim == T ]; then
  rm -r "$simFold"parameter_sensitivity_files
fi


# end cleanup #
