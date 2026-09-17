#!/bin/bash

#SBATCH --hint=nomultithread
#SBATCH --job-name=ecland_IAV
#SBATCH --output=log/ecland_IAV-%j.out
#SBATCH --error=log/ecland_IAV-%j.err
#SBATCH --time=2:00:00

# ============================================
#    Script create-time-varying-veg.sh
#    -----------------------------------
#
#    DESCRIPTION:
#    Create time varying LC and LAI
#
#    USAGE :
#    create-time-varying-veg.sh [GROUP] [REGION] [STARTYEAR] [ENDYEAR] [CLIM]
#
#    EXAMPLES:
#    sbatch create-time-varying-veg.sh YEARLY_2D iberia 2015 2019 true
#
#    ARGUMENTS:
#      * GROUP: Group name related to the lenght of the run:
#             - YEARLY_2D
#             - MONTHLY_2D
#             - REG25
#             - REG10
#      * REGION: Region of the Earth:
#             - iberia
#             - amazon
#             - europe
#             - southamerica
#      * STARTYEAR: Initial year of the simulation: from 1999 to 2019
#      * ENDYEAR: Last year of the simulation: from 1999 to 2019
#      * CLIM: true (LAI 2010-2019 clim CGLS) or false (LAI year CGLS)
#    Author:
#    --------
#    I. Ayan-Miguez (BSC) January 2026. Based on the surface_model.ksh script from ECMWF.
#
#=============================================

set -ue

#Default
export GROUP=$1
export REGION=$2
export STARTYEAR=$3
ENDYEAR=$4
CLIM=${5:-false}

source ${PATH_source}/scripts/env.sh

set -xv

#Define directories
SCRIPTDIR=${PATH_source}/scripts/time-varying-LC-LAI
RUNDIR=${ECLAND_TMP}/time-varying-veg/${REGION}_${STARTYEAR}_${ENDYEAR}
export INIDATA=${ECLAND_REPO}/era5/clim

YYYY=$STARTYEAR

rm -rf $RUNDIR ; mkdir -p $RUNDIR
cp -f ${SCRIPTDIR}/*.* $RUNDIR
cd ${RUNDIR}

ncks -O -v landsea ${INIDATA}/${GROUP}/surfclim_${REGION}_${STARTYEAR}01-* lsm_tmp.nc
ncap2 -O -s 'lon=lon-360*(lon>180)' lsm_tmp.nc lsm_regular.nc

rm *tmp.nc

export HDF5_USE_FILE_LOCKING=FALSE

cat > setup_cython.py << EOF
from distutils.core import setup
from Cython.Build import cythonize
import numpy
setup(
  name = 'cython ext',
  ext_modules = cythonize("land_cover_veg_cython_utils.pyx"),
  include_dirs=[numpy.get_include(),],
  )
EOF

rm -rf build land_cover_veg_cython_utils.c  land_cover_veg_cython_utils.cpython-38-x86_64-linux-gnu.so # for security to be sure the cython extension is always re-build

python3 setup_cython.py build_ext --inplace

while [[ $YYYY -le $ENDYEAR ]] ; do

  if [[ $CLIM == true ]]; then
    rm -rf clim_${YYYY}; mkdir -p clim_${YYYY}
    ./update_IAV_veg.sh -n clim_${YYYY} -t ${YYYY} -g lsm_regular.nc
    cp clim_${YYYY}/surfviav_${YYYY} ${INIDATA}/${GROUP}/surfviav_${REGION}_clim_${YYYY}.nc
  else
    rm -rf ${YYYY}; mkdir -p ${YYYY}
    ./update_IAV_veg.sh -n ${YYYY} -t ${YYYY} -g lsm_regular.nc
    cp ${YYYY}/surfviav_${YYYY} ${INIDATA}/${GROUP}/surfviav_${REGION}_${YYYY}.nc
  fi

  YYYY=`expr $YYYY + 1`

done

rm -rf ${RUNDIR}