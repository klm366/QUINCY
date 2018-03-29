#!/bin/bash

# this script:
# runs at the end of each day at 11 pm via cronjob (crontab) at io1 (Tatra (IO1/Galactica))
# uses an extra local working copy of SLM/trunk at scratch/jengel/autoslmtestrun/
# runs with jengel user rights
#

# global variables
autoTestRunFolder="autoslmtestrun/"
bgiJEScratchFolder="/Net/Groups/BGI/scratch/jengel/"
workingCopyFolder="$bgiJEScratchFolder""$autoTestRunFolder"trunk/


###################################################################################################
## after all the cleanup-jobs have ended the postprocessing job may execute this code            ##
###################################################################################################

if [ "$1" == "postprocessing" ]
  then

    svnRev="$2"
    jobNameD="$3"
    jobNameD2="$4"
    jobNameS="$5"
    jobNameS2="$6"
    jobNameS3="$7"
    jobNameS4="$8"

    echo " starting postprocessing of auto SLM testruns for Revision: ""$svnRev"
    date
    echo " "
    echo " jobs are:"
    echo "  ""$jobNameD"
    echo "  ""$jobNameD2"
    echo "  ""$jobNameS"
    echo "  ""$jobNameS2"
    echo "  ""$jobNameS3"
    echo "  ""$jobNameS4"
    echo " "

    saveFolder="$bgiJEScratchFolder""$autoTestRunFolder""$svnRev""/"
    mkdir -p "$saveFolder"
    sync

    ## do some stats and analysis on the debug-run simulation results ##
    echo "analysis of: ""$jobNameD"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameD" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameD"/*/error/err* > "$saveFolder""$jobNameD"_err.txt
    sync
    ls -lah "$saveFolder""$jobNameD"_err.txt
    echo " "

    ## do some stats and analysis on the debug-run simulation results ##
    echo "analysis of: ""$jobNameD2"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameD2" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameD2"/*/error/err* > "$saveFolder""$jobNameD2"_err.txt
    sync
    ls -lah "$saveFolder""$jobNameD2"_err.txt
    echo " "

    ## do some stats and analysis on the standard-run simulation results ##
    echo "analysis of: ""$jobNameS"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameS" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameS"/*/error/err* > "$saveFolder""$jobNameS"_err.txt
    sync
    ls -lah "$saveFolder""$jobNameS"_err.txt
    echo " "

    ## do some stats and analysis on the standard-run simulation results ##
    echo "analysis of: ""$jobNameS2"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameS2" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameS2"/*/error/err* > "$saveFolder""$jobNameS2"_err.txt    
    sync
    ls -lah "$saveFolder""$jobNameS2"_err.txt
    echo " "

    ## do some stats and analysis on the standard-run simulation results ##
    echo "analysis of: ""$jobNameS3"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameS3" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameS3"/*/error/err* > "$saveFolder""$jobNameS3"_err.txt
    sync
    ls -lah "$saveFolder""$jobNameS3"_err.txt
    echo " "

    ## do some stats and analysis on the standard-run simulation results ##
    echo "analysis of: ""$jobNameS4"
    # move simulation results
    mv "$bgiJEScratchFolder"slm/"$jobNameS4" "$saveFolder"
    sync && sleep 2
    # extract potential error messages into one file
    cat "$saveFolder""$jobNameS4"/*/error/err* > "$saveFolder""$jobNameS4"_err.txt    
    sync
    ls -lah "$saveFolder""$jobNameS4"_err.txt
    echo " "    

  # ok postprocessing done exiting this script
  echo " postprocessing of automated testruns finished"
  date
  echo " "
  echo " "
  exit 0
fi


###################################################################################################
## start the simulations & cleanup jobs & postprocessing job                                     ##
###################################################################################################


# cd into local working copy
cd "$workingCopyFolder"

# revision number of repository trunk/src/ and local working copy
revRepoTrunk=$(svn info https://projects.bgc-jena.mpg.de/QUINCY/svn/project-A/trunk | grep "Rev:" | awk '{print $4}')
localWCopy=$(svn info ./ | grep "Rev:" | awk '{print $4}')
revRepoTrunkSrc=$(svn info https://projects.bgc-jena.mpg.de/QUINCY/svn/project-A/trunk/src | grep "Rev:" | awk '{print $4}')
localWCopySrc=$(svn info ./src | grep "Rev:" | awk '{print $4}')

# run test runs only if a newer revision of the trunk/src/ folder is available | or if "ignoreRev" is passed as argument to the script (for manual testing purpose)
if [ $revRepoTrunkSrc -gt $localWCopySrc ] || [ "$1" == "ignoreRev" ]
  then
  
  # check for changing the customSiteList= F/T, per default it is T
  # check whether the "allsites" argument is present, then it will change to F
  flagCustomSiteList="T"
  allSitesArg="$(grep -w "allsites" <<< "$@")"  # check the list of arguments for "allsites" string
  if [ ! -z "${allSitesArg// }" ]               # -z test is explicitly for: the length of string is zero
   then 
    flagCustomSiteList="F"
    echo " autotestrun simulations will be done for all fluxnet sites, i.e., not using the auto_testrun_flux_sites_list.dat."
    echo " "
  fi    
  
	# update local working copy
	svn update && sleep 2

  ##########################
  # compile in debug mode
  cd "$workingCopyFolder"build/
  rm -f Make.inc.local && ln -s Make.inc.local.ifort_debug_BGC Make.inc.local && make clean && sync && make -j 4
  sleep 2

  cd ../contrib/site_level/
  
  # run 100 years fluxnet simulations in debug mode with simple_0d: usecase 02
  jobNameD="auto_r""$revRepoTrunk""_100yfluxdef_dbg_s0d"
  ./prep_run_slm.sh -u 02 -e "$jobNameD" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameD="exp_""$jobNameD"
  sleep 30s 
  # run 100 years fluxnet simulations in debug mode with simple_1d: usecase 04
  jobNameD2="auto_r""$revRepoTrunk""_100yfluxdef_dbg_s1d"
  ./prep_run_slm.sh -u 04 -e "$jobNameD2" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameD2="exp_""$jobNameD2"

  ##########################
  # compile in performance mode
  cd "$workingCopyFolder"build/
  rm -f Make.inc.local && ln -s Make.inc.local.ifort_perform_BGC Make.inc.local && make clean && sync && make -j 4
  sleep 2

  cd ../contrib/site_level/
  
  # run 500 years fluxnet simulations in standard mode with simple_0d: usecase 03
  jobNameS="auto_r""$revRepoTrunk""_500yfluxdef_std_s0d"
  ./prep_run_slm.sh -u 03 -e "$jobNameS" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameS="exp_""$jobNameS"
  sleep 30s 
  # run 500 years fluxnet simulations in standard mode with simple_1d: usecase 05
  jobNameS2="auto_r""$revRepoTrunk""_500yfluxdef_std_s1d"
  ./prep_run_slm.sh -u 05 -e "$jobNameS2" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameS2="exp_""$jobNameS2"
  sleep 30s 
  # run 500 years fluxnet simulations in standard mode with SLM Soil jsm: usecase 06
  jobNameS3="auto_r""$revRepoTrunk""_500yfluxdef_std_soiljsm"
  ./prep_run_slm.sh -u 06 -e "$jobNameS3" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameS3="exp_""$jobNameS3"
  sleep 30s 
  # run 500 years fluxnet simulations in standard mode with simple_1d and P cycle: usecase 07
  jobNameS4="auto_r""$revRepoTrunk""_500yfluxdef_std_s1d_p_cycle"
  ./prep_run_slm.sh -u 07 -e "$jobNameS4" -c "$flagCustomSiteList" -n auto_testrun_flux_sites_list.dat
  jobNameS4="exp_""$jobNameS4"


  ####################
  # jobs are running #
  ####################


  ##########################
  # start the postprocessing job that calls this script and executes the above postprocessing code only
  sleep 30
  echo " "
  echo " preparing and sending to the queue the postprocessing job"
      
  jobD="$jobNameD""Cleanup"
  jobD2="$jobNameD2""Cleanup"
  jobS="$jobNameS""Cleanup"
  jobS2="$jobNameS2""Cleanup"
  jobS3="$jobNameS3""Cleanup"
  jobS4="$jobNameS4""Cleanup"
  jobsToWait=$(echo "ended(""$jobD"") && ended(""$jobD2"") && ended(""$jobS"") && ended(""$jobS2"") && ended(""$jobS3"") && ended(""$jobS4"")")
  bsub -w "$jobsToWait" -q "io" -J "autotestrun_post" "$workingCopyFolder""contrib/additional_scripts/automated_slm_testrun.sh" postprocessing "$revRepoTrunk" "$jobNameD" "$jobNameD2" "$jobNameS" "$jobNameS2" "$jobNameS3" "$jobNameS4" 

  echo " "
  echo "Done!"

  #############################################################################
  #############################################################################
else  # if [ $revRepoTrunk -gt $localWCopy ]
    cd "$bgiJEScratchFolder""$autoTestRunFolder"
    # check if a log file already exists, if not create one (sed relies on a lof file beeing present)
    if [ ! -e ./cronjob-slm-autotestrun.log ]
     then
      echo " " > ./cronjob-slm-autotestrun.log
      sync
    fi
    dateNow=$(date "+%d.%m.%Y %X")
    new=$(echo "$dateNow"" | repo ./src : ""$revRepoTrunkSrc"" | local ./src : ""$localWCopySrc")
    sed -i 1i"$new" ./cronjob-slm-autotestrun.log    # inserts $new in line number 1
    sed -i '101,110d' ./cronjob-slm-autotestrun.log  # deletes all lines 101-110 (100 lines remain); no error if there are no lines 102-110
fi


# fin #
