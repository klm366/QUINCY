#!/bin/sh

cd $T/git/jsbach/jsbach-dev/src
rsync -Rav `git status -s | grep '^ M ' | awk '{print $2}'` thunder4:/scratch/mpi/mj0060/m212005/icon-dev-postswitch/src/lnd_phy_jsbach/

