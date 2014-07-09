#!/bin/ksh
#------------------------------------------------------------------------------
# Script to construct from a given ICON grid file
# initial fields and boundary conditions consistent
# for the land model JSBACH4+ and the ICON atmosphere.
# 1. compile and run jsbach_init_file.f90 via the script jsbach_init_file.ksh to construct
# a JSBACH input file on a gaussian grid, which can also be used for ECHAM6
# 2. run gauss_to_icon.ksh to interpolate the fields to the ICON grid
#
# jsbach_init_file.f90 was coded by Thomas Raddatz and Veronika Gayler at MPI-Met
# gauss_to_icon.ksh combines a script coded by Veronika Gayler (for land fields) and one
# coded by Marco Giorgetta (for atmosphere fields).
#
# Thomas Raddatz 4. April 2013

set -e

set +u
. ${MODULESHOME}/init/ksh
set -u
# control fortran compiler (use version 5.3.907)
module unload nag
module load nag/5.3.907
# control cdo version (use version 1.6.0)
module unload cdo
module load cdo/1.6.0

# get the script jsbach_init_file.ksh and the program jsbach_init_file.f90
#svn checkout -r 6285 https://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/contrib/initial_tarfiles
svn checkout https://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/contrib/initial_tarfiles
mv ./initial_tarfiles/jsbach_init_file.ksh .
mv ./initial_tarfiles/jsbach_init_file.f90 .

# switch off interactive mode in the script jsbach_init_file.ksh
export interactive_ctrl=1

# component model resolutions
export res_atm=T63         # atmosphere horizontal resolution
export res_oce=GR15        # ocean horizontal resolution

# ICON grid file
pool_grid=/pool/data/ICON/grids/private/r2b4_amip # pool with ICON grid data
export grid=r2b4_amip
export gridfile=${pool_grid}/${grid}.nc

# echam6 input data
export pool=/pool/data/ECHAM6/input/r0003/${res_atm}  # pool with atm input data
export pool_land=/pool/data/JSBACH/prepare/${res_atm} # pool with srf input data

# script directory
startdir=$(pwd)
[[ -f gauss_to_icon.ksh ]] || cp `dirname ${0}`/gauss_to_icon.ksh .

# parameters for jsbach_init_file.f90
export dynveg=false
export year_ct=1976; export year_cf=1976
export ntiles=11
export landcover_series=false
export c3c4crop=true
export cmip5_pasture=true
export l5_soil_layers=true
export lroot_depth=${l5_soil_layers}
export echam_fractional=false
export srcdir=${startdir}
export wrkdir=${startdir}

# atmosphere input
mkdir -p input/echam6
cd input/echam6
cp -p ${pool}/${res_atm}${res_oce}_VGRATCLIM.nc  .
cp -p ${pool}/${res_atm}${res_oce}_VLTCLIM.nc    .
cp -p ${pool}/${res_atm}${res_oce}_jan_surf.nc   .
cp -p ${pool}/${res_atm}_TSLCLIM2.nc             .
cd ../..

# land (JSBACH) input
mkdir -p ${srcdir}/input/jsbach/${grid}
cd ${srcdir}/input/jsbach/${grid}

${srcdir}/jsbach_init_file.ksh

# remove the executable and modules
rm jsbach_init_file mo_kinds.mod mo_vegparams.mod

${srcdir}/gauss_to_icon.ksh

exit
