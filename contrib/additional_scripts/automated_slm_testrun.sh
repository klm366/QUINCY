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
###################################################################################################

# cd into local working copy
cd "$workingCopyFolder"

# revision number of repository trunk/src/ and local working copy
revRepoTrunk=$(svn info https://projects.bgc-jena.mpg.de/QUINCY/svn/project-A/trunk/src | grep "Rev:" | awk '{print $4}')
localWCopy=$(svn info ./src | grep "Rev:" | awk '{print $4}')

# run test runs only if a newer revision of the trunk/src/ folder is available
if [ $revRepoTrunk -gt $localWCopy ]
  then
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
  ./prep_run_slm.sh -u 02 -e "$jobNameD"
  jobNameD="exp_""$jobNameD"
  # run 100 years fluxnet simulations in debug mode with simple_1d: usecase 04
  jobNameD2="auto_r""$revRepoTrunk""_100yfluxdef_dbg_s1d"
  ./prep_run_slm.sh -u 04 -e "$jobNameD2"
  jobNameD2="exp_""$jobNameD2"

  ##########################
  # compile in performance mode
  cd "$workingCopyFolder"build/
  rm -f Make.inc.local && ln -s Make.inc.local.ifort_perform_BGC Make.inc.local && make clean && sync && make -j 4
  sleep 2

  cd ../contrib/site_level/
  
  # run 300 years fluxnet simulations in standard mode with simple_0d: usecase 03
  jobNameS="auto_r""$revRepoTrunk""_300yfluxdef_std_s0d"
  ./prep_run_slm.sh -u 03 -e "$jobNameS"
  jobNameS="exp_""$jobNameS"
  # run 300 years fluxnet simulations in standard mode with simple_1d: usecase 05
  jobNameS2="auto_r""$revRepoTrunk""_300yfluxdef_std_s1d"
  ./prep_run_slm.sh -u 05 -e "$jobNameS2"
  jobNameS2="exp_""$jobNameS2"


  ####################################
  # get the job ID of the cleanup jobs
  sleep 30
  jobIDCleanD=$(bjobs -J "$jobNameD""Cleanup" | tail -n 1 | awk '{print $1}')
  jobIDCleanD2=$(bjobs -J "$jobNameD2""Cleanup" | tail -n 1 | awk '{print $1}')
  jobIDCleanS=$(bjobs -J "$jobNameS""Cleanup" | tail -n 1 | awk '{print $1}')
  jobIDCleanS2=$(bjobs -J "$jobNameS2""Cleanup" | tail -n 1 | awk '{print $1}')

  # return to working copy folder
  cd "$workingCopyFolder"


  ####################
  # jobs are running #
  ####################

  ## debug - usecase 02 #######################################################
  # stay in the while loop until the cleanup-job of the debug-run has finished
  status="RUN"
  countExit=1
  while [ $status == "RUN" ] || [ $status == "PEND" ]
    do
      sleep 3m
      status=$(bjobs "$jobIDCleanD" | tail -n 1 | awk '{print $3}')
      ((countExit++))
      # change queue to io if 5 hours have passed
      if [ $countExit -eq 100 ]; then bswitch io "$jobIDCleanD"; fi
      # exit this script if more than 10 hours have passed
      if [ $countExit -gt 200 ]
       then
        echo "The status of jobIDCleanD is ""$status"" and countExit is ""$countExit"
        date "+%d.%m.%Y %X"
        echo "EXIT !"
        exit 1
      fi
  done
  ## do some stats and analysis on the debug-run simulation results ##
  date
  echo "do some stats and analysis on the debug-run simulation results"
  echo " "
  # move simulation results
  mv "$bgiJEScratchFolder"slm/"$jobNameD" "$bgiJEScratchFolder""$autoTestRunFolder"
  sync && sleep 2
  # extract potential error messages into one file
  cat "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameD"/*/error/err* > "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameD"_err.txt

  ## debug - usecase 04 #######################################################
  # stay in the while loop until the cleanup-job of the debug-run has finished
  #  but first check whether the 'bjobs' command does finish with exit 0
  #   (note: in case the job (simulations) already finished a while ago, the bjobs command would not show any status information
  #    about this job anymore (like RUN | PEND) but return with an error message that there is no info about this job)
  sleep 10
  bjobs "$jobIDCleanD2" > /dev/null
  returnB=$(echo $?)
  if [ "$returnB" -eq 0 ] # if bjobs does return error-free
   then
    status="RUN"
    countExit=1           # if bjobs returns an error code (i.e. job long finished), set 'status' to any but "RUN" or "PEND" to not enter the while loop
   else
    status="JOB_ALREADY_FINISHED"
  fi
  # if job is still running/peding enter the while loop, waiting for it to finish
  while [ $status == "RUN" ] || [ $status == "PEND" ]
    do
      sleep 3m
      status=$(bjobs "$jobIDCleanD2" | tail -n 1 | awk '{print $3}')
      ((countExit++))
      # change queue to io if additional 2 1/2 hours have passed
      if [ $countExit -eq 50 ]; then bswitch io "$jobIDCleanD2"; fi
      # exit this script if additional 5 hours have passed
      if [ $countExit -gt 100 ]
       then
        echo "The status of jobIDCleanD2 is ""$status"" and countExit is ""$countExit"
        date "+%d.%m.%Y %X"
        echo "EXIT !"
        exit 1
      fi
  done
  ## do some stats and analysis on the debug-run simulation results ##
  date
  echo "do some stats and analysis on the debug-run simulation results"
  echo " "
  # move simulation results
  mv "$bgiJEScratchFolder"slm/"$jobNameD2" "$bgiJEScratchFolder""$autoTestRunFolder"
  sync && sleep 2
  # extract potential error messages into one file
  cat "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameD2"/*/error/err* > "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameD2"_err.txt

  ## standard - usecase 03 #######################################################
  # stay in the while loop until the cleanup-job of the standard-run has finished
  sleep 10
  bjobs "$jobIDCleanS" > /dev/null
  returnB=$(echo $?)
  if [ "$returnB" -eq 0 ] # if bjobs does return error-free
   then
    status="RUN"
    countExit=1           # if bjobs returns an error code (i.e. job long finished), set 'status' to any but "RUN" or "PEND" to not enter the while loop
   else
    status="JOB_ALREADY_FINISHED"
  fi
  # if job is still running/peding enter the while loop, waiting for it to finish
  while [ $status == "RUN" ] || [ $status == "PEND" ]
    do
      sleep 3m
      status=$(bjobs "$jobIDCleanS" | tail -n 1 | awk '{print $3}')
      ((countExit++))
      # change queue to io if additional 2 1/2 hours have passed
      if [ $countExit -eq 50 ]; then bswitch io "$jobIDCleanS"; fi
      # exit this script if additional 5 hours have passed
      if [ $countExit -gt 100 ]
       then
        echo "The status of jobIDCleanS is ""$status"" and countExit is ""$countExit"
        date "+%d.%m.%Y %X"
        echo "EXIT !"
        exit 1
      fi
  done
  ## do some stats and analysis on the standard-run simulation results ##
  date
  echo "do some stats and analysis on the standard-run simulation results"
  echo " "
  # move simulation results
  mv "$bgiJEScratchFolder"slm/"$jobNameS" "$bgiJEScratchFolder""$autoTestRunFolder"
  sync && sleep 2
  # extract potential error messages into one file
  cat "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameS"/*/error/err* > "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameS"_err.txt

  ## standard - usecase 05 #######################################################
  # stay in the while loop until the cleanup-job of the standard-run has finished
  sleep 10
  bjobs "$jobIDCleanS2" > /dev/null
  returnB=$(echo $?)
  if [ "$returnB" -eq 0 ] # if bjobs does return error-free
   then
    status="RUN"
    countExit=1           # if bjobs returns an error code (i.e. job long finished), set 'status' to any but "RUN" or "PEND" to not enter the while loop
   else
    status="JOB_ALREADY_FINISHED"
  fi
  # if job is still running/peding enter the while loop, waiting for it to finish
  while [ $status == "RUN" ] || [ $status == "PEND" ]
    do
      sleep 3m
      status=$(bjobs "$jobIDCleanS2" | tail -n 1 | awk '{print $3}')
      ((countExit++))
      # change queue to io if additional 2 1/2 hours have passed
      if [ $countExit -eq 50 ]; then bswitch io "$jobIDCleanS2"; fi
      # exit this script if additional 5 hours have passed
      if [ $countExit -gt 100 ]
       then
        echo "The status of jobIDCleanS2 is ""$status"" and countExit is ""$countExit"
        date "+%d.%m.%Y %X"
        echo "EXIT !"
        exit 1
      fi
  done
  ## do some stats and analysis on the standard-run simulation results ##
  date
  echo "do some stats and analysis on the standard-run simulation results"
  echo " "
  # move simulation results
  mv "$bgiJEScratchFolder"slm/"$jobNameS2" "$bgiJEScratchFolder""$autoTestRunFolder"
  sync && sleep 2
  # extract potential error messages into one file
  cat "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameS2"/*/error/err* > "$bgiJEScratchFolder""$autoTestRunFolder""$jobNameS2"_err.txt


  #############################################################################
  else
    cd "$bgiJEScratchFolder""$autoTestRunFolder"
    # check if a log file already exists, if not create one (sed relies on a lof file beeing present)
    if [ ! -e ./cronjob-slm-autotestrun.log ]
     then
      echo " " > ./cronjob-slm-autotestrun.log
      sync
    fi
    dateNow=$(date "+%d.%m.%Y %X")
    new=$(echo "$dateNow"" | repo: ""$revRepoTrunk"" | local: ""$localWCopy")
    sed -i 1i"$new" ./cronjob-slm-autotestrun.log    # inserts $new in line number 1
    sed -i '101,110d' ./cronjob-slm-autotestrun.log  # deletes all lines 101-110 (100 lines remain); no error if there are no lines 102-110
fi


# fin #

