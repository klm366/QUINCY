#!/bin/bash

## 
##
## script to cleanup temporary files and folders after a whole job-array of slm simulations was completed
## called from the prep_run_slm.sh
## this script will not be copied anywhere but remain at its place (contrib/site_level) and work from here


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################

## arguments passed to cleanup_slm_hpc_run.sh ##
jobName=$1                        # a user specific alphanumeric string
array=$2                          # was an job array or just one job simulated  T/F
# $3 == singleSiteName provided only if just one single site was run


## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"
slmModel=$(sed -n 3p "$simFold""default_files/config.txt")                                    # SLM model: land plant canopy soil
paramSensitivitySim=$(sed -n 7p "$simFold""default_files/config.txt" | awk '{print $1}')      # T/F  -> parameter sensitivity
compareSimCRUNCEPvsGFDB=$(sed -n 8p "$simFold""default_files/config.txt" | awk '{print $1}')  # T/F  -> compare with global forest data base


## move job specific error file into the particular job folder ##
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
          mv "$simFold""error/err.""$errFileCount"".txt"  "$simFold""$simCount"_*/error
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
    # move error files
    for i in "$simFold"error/*
     do
      siteNumber=$(awk -F "." '{print $2}' <<< $i)
      siteName=$(sed -n "$siteNumber"p "$siteList" | awk '{print $1}')
      mv "$simFold""error/err.""$siteNumber"".txt" "$simFold""$siteName""/error/"
    done
  fi
 # $array == no
 else
  singleSiteName=$3               # name of the site
  mv "$simFold""error/err.txt" "$simFold""$singleSiteName""/error/"
fi

## create an md5sum file for all output files, i.e. files in the folder "output"
md5sum "$simFold"*/output/* > "$simFold"md5sum-output-files.txt

## save the binary and the parameter table
mkdir -p "$simFold"x_binary
mv "$simFold"default_files/*.x "$simFold"x_binary
if [ $paramSensitivitySim == T ]; then
  mkdir -p "$simFold"parameter_list
  mv "$simFold"parameter_sensitivity_files/param_table.txt "$simFold"parameter_list
fi

## run two R Scripts for comparing simulated CRUNCEP data with global forest data base data at daily basis
if [ $compareSimCRUNCEPvsGFDB == T ]; then
  Rscript "$simFold"default_files/sim_cruncep_vs_gfdb-process_daily.R $simFold  # creates two .csv files
  sync
  Rscript "$simFold"default_files/sim_cruncep_vs_gfdb-create_figs.R $simFold    # creates a .pdf with multiple figures
  rm "$simFold"obs_site.csv
  sync
fi

## remove some folders not needed anymore ##
rm -r "$simFold"default_files
rm -r "$simFold"error
rm -r "$simFold"out
if [ $paramSensitivitySim == T ]; then
  rm -r "$simFold"parameter_sensitivity_files
fi

## create one pdf file from the postprocessing pdf of all sites
if [ "$slmModel" == "land" ] || [ "$slmModel" == "plant" ] # || [ "$slmModel" == "soil" ]
 then

  cd "$simFold"
  mkdir -p ./pdftmp
  count=100

  for site in $(ls -d ./*/)
   do

    siteBase=$(basename $site)  # just: HAR_01 SAS_04 TAP_67

    if [ $siteBase != "x_binary" ] && [ $siteBase != "pdftmp" ]
    then
      ((count++))

      cp ./"$siteBase"/postprocessing/*.pdf  ./pdftmp/"$count"b.pdf

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
fi

# end cleanup #
