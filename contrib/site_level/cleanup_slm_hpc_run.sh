#!/bin/bash

## 20160601
##
## script to cleanup temporary files and folders after a whole job-array of slm simulations was completed
## called from the prep_run_slm.sh
## this script will not be copied anywhere but remain at its place (contrib/site_level) and work from here


########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################

## arguments passed to cleanup_slm_hpc_run.sh ##
jobName=$1                        # a user specific alphanumeric string
array=$2                          # was an job array or just one job simulated


## basic folders ##
username=$(whoami)
bgiScratchUserFold="/Net/Groups/BGI/scratch/""$username""/slm/"
simFold="$bgiScratchUserFold""$jobName""/"


## move job specific error file into the particular job folder ##
if [ "$array" == "yes" ]
 then
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
 else # i.e. $array == no
  singleSiteName=$3               # name of the site
  mv "$simFold""error/err.txt" "$simFold""$singleSiteName""/error/"
fi

## create an md5sum file for all output files, i.e. files in the folder "output"
md5sum "$simFold"*/output/* > "$simFold"md5sum-output-files.txt


## remove some folders not needed anymore ##
rm -r "$simFold"default_files
rm -r "$simFold"error
rm -r "$simFold"out



