#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2026, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=aggregate_LUH_transitions.o%j
#SBATCH --error=aggregate_LUH_transitions.o%j
#SBATCH --account=bm0891
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=5120
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --ntasks=1

#----------------------------------------------------------------------
# AKN April 2026, script to aggregate LUH type transitions
# - With rangelands treated depending on the fstnf value (i.e., potential vegetation type) from the static file
#   - if fstnf = 1, i.e. forests in LUH, transitions to/from rangelands are treated as transtions to/from pasture;
#   - if fstnf = 0, i.e., no forest in LUH, transitions to/from rangelands are treated as transitions to/from natural vegetation;
# - Transitions from/to urban are treated as transition from/to pasture. However, for the states 'urban' is assigned to bare land.
# - Negative transitions rates are set to zero. They seem to occur occationally for transition rates to/from c3nfx.


start_year=1979
end_year=1982
usecase="12pfts"  # 12pfts, 11pfts

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

scriptName=aggregate_LUH_transitions.sh

# usually the original grid (for further use in the icon-land init file scripts)
targetGrid="global_0.25"
# -> note: lats in LUH data have another orientation as compared to the one in global_0.25
#          therefore a remapping here makes sense to ensure that all 0.25 grids used have the same orientation!

tag="CMIP7_forcings"
version="v3.1.1"

mainTransFile=/work/bb1469/b383102/${tag}/RAW_DATA/land/${version}/multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2023.nc
mainStatesFile=/work/bb1469/b383102/${tag}/RAW_DATA/land/${version}/multiple-states_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2024.nc
fixedInputFile=/work/bb1469/b383102/${tag}/RAW_DATA/land/${version}/multiple-static_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn.nc

#outfolder=./YOUROUTPUTFOLDER
outfolder=./test_output

#-------------------------------------
echo "--> ${scriptName}"

# Some global attributes to be set
luh3_preprocessing="Aggregated LUH transitions (using ${scriptName}) and ${tag} LUH data ${version}"
luh3_preprocessing="${luh3_preprocessing}; contact: a.naumann@lmu.de"
luh3_comment="Urban is treated as pasture."
luh3_comment="${luh3_comment} Rangelands on potential forests are treated as pasture (LUH static data: fstnf = 1), else: natural vegetation"

# mapping of original to aggregated states
# NOTE: Range is listed twice. Rangelands are aggregated to pasture or natural vegetation depending on the fstnf value
org_states=( primf primn urban secdf secdn pastr range range c3ann c3per c4ann c4per c3nfx )
if [ "$usecase" = "11pfts" ]; then
  aggregated_states=( natural natural pasture natural natural pasture natural pasture crop crop crop crop crop )
elif [ "$usecase" = "12pfts" ]; then
  aggregated_states=( natural natural pasture natural natural pasture natural pasture c3crop c3crop c4crop c4crop c3crop )
else
  echo "unknown usecase: ${usecase}"
  exit 1
fi
# Note: the "pltns" state variable and accoridng transitions contain zeros only in LUH3 V1; prop. needs to be handled > V1

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

# get fstnf file
fstnfFile=./multiple-static_input4MIPs_fstnf.nc
${cdo} -selvar,fstnf ${fixedInputFile} ${fstnfFile}

# get selected years from multiple states and transitions file
statesFile=./multiple-states_input4MIPs_landState_${start_year}-${end_year}.nc
${cdo} -selyear,${start_year}/${end_year} ${mainStatesFile} ${statesFile}
transFile=./multiple-transitions_input4MIPs_landState_${start_year}-${end_year}.nc
${cdo} -selyear,${start_year}/${end_year} ${mainTransFile} ${transFile}


#----------------------
#-- aggregate transitions
nonExistingTransitions=( primf_to_secdf primn_to_secdn range_to_range )
# and there are obviously no transitions to primf or primn

state_from_index=0
fileList=""
while [ ${state_from_index} -lt ${#org_states[@]} ]; do
  state_to_index=0
  state_from=${org_states[$state_from_index]}
  aggregated_state_from=${aggregated_states[$state_from_index]}
  while [ ${state_to_index} -lt ${#org_states[@]} ]; do
    state_to=${org_states[$state_to_index]}
    aggregated_state_to=${aggregated_states[$state_to_index]}

    thisVar="${state_from}_to_${state_to}"
    echo "current var: ${thisVar}"

    if [ ! ${aggregated_state_from} = ${aggregated_state_to} ]; then
      thisTargetVar="${aggregated_state_from}2${aggregated_state_to}"

      # there are no transitions to primf or primn
      if [ ${state_to} = primf ] || [ ${state_to} = primn ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      doesNotExist=false
      for nonTrans in ${nonExistingTransitions[@]}; do
        if [ ${nonTrans} = ${thisVar} ]; then
          doesNotExist=true
          break
        fi
      done
      if [ $doesNotExist = true ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      echo "... adding ${thisVar} to ${thisTargetVar}"

      # get the data
      ${cdo} -selvar,${thisVar} ${transFile} tmp1.nc

      # treat rangelands
      if [ ${state_from} = range ] && [ ${aggregated_state_from} = natural ]; then
          # only take those transitions where LUH does not assume forest (fstnf = 0)
          ${cdo} -mul tmp1.nc -mulc,-1 -subc,1 ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_to} = range ] && [ ${aggregated_state_to} = natural ]; then
          ${cdo} -mul tmp1.nc -mulc,-1 -subc,1 ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_to} = range ] && [ ${aggregated_state_to} = pasture ]; then
          # only take those transitions where LUH assumes forest (fstnf = 1)
          ${cdo} -mul tmp1.nc ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_from} = range ] && [ ${aggregated_state_from} = pasture ]; then
          ${cdo} -mul tmp1.nc ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      fi

      # set values <0 to 0, as there are some negative values in the original data
      ${cdo} -expr,"${thisVar}=(${thisVar}<0)?0:${thisVar}" tmp1.nc tmp2.nc
      mv tmp2.nc tmp1.nc

      thisAggregationFile=tmp_${thisTargetVar}.nc
      if [ -f ${thisAggregationFile} ]; then
        ${cdo} -add ${thisAggregationFile} tmp1.nc tmp2.nc
        mv tmp2.nc ${thisAggregationFile}
      else
        ${cdo} -chname,${thisVar},${thisTargetVar} tmp1.nc ${thisAggregationFile}
        # remove old not suitable attributes
        ${ncatted} -a standard_name,${thisTargetVar},d,, ${thisAggregationFile}
        fileList="${fileList} ${thisAggregationFile}"
      fi

      rm -f tmp.nc tmp1.nc tmp2.nc
    fi
    (( state_to_index = state_to_index + 1 ))
  done

  (( state_from_index = state_from_index + 1 ))
done

#-------------------
# join the files and remove lat/lon bounds
${cdo} merge ${fileList} tmp_all_transitions.nc
rm ${fileList}
ncks -C -O -x -v lat_bnds,lon_bnds tmp_all_transitions.nc tmp_all_transitions.nc

# remap to target grid
${cdo} -remapcon,${targetGrid} tmp_all_transitions.nc LUH_trans_${start_year}-${end_year}.nc
rm tmp_all_transitions.nc

# set global attributes
${ncatted} -O -a ,global,d,, LUH_trans_${start_year}-${end_year}.nc
${ncatted} -O -a luh3_comment,global,o,c,"${luh3_comment}" LUH_trans_${start_year}-${end_year}.nc
${ncatted} -O -a luh3_preprocessing,global,o,c,"${luh3_preprocessing}" LUH_trans_${start_year}-${end_year}.nc


#------------------
# remapping to ICON grid and rescaling
# - In LUH3, land-use transitions are given as the fraction of each grid-cell that transitions one from land-use
#   state to another in a given year. JSBACH needs the transition rates relative to the vegetated fraction of the grid cell.
#   Therefore, the transition rates from LUH3 need to be scaled with land * fract_veg * veg_ratio_max.
# --- TBD: Do remapping to ICON grid and rescaling here only for test purposes.
#          Later this should be incorporated in the ICON-Land preprocessing chain, possibly in generate_bc_land_fraction.sh.

testing_only=false
if [ "$testing_only" = true ]; then
  echo "Remapping to ICON grid for testing purpose..."
  targetGrid2=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc   # mpim R2B4 grid
  # targetGrid2=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc   # dwd R2B5 grid
  ${cdo} -remapcon,${targetGrid2} LUH_trans_${start_year}-${end_year}.nc tmp_remapped.nc
  rm LUH_trans_${start_year}-${end_year}.nc

  echo "Scaling transitions with land * fract_veg * veg_ratio_max..."
  # vrm is constant in time
  if [ "$usecase" = "11pfts" ]; then
    vrm_file=/pool/data/ICON/grids/public/mpim/0049/land/r0101/bc_land_frac.nc
  elif [ "$usecase" = "12pfts" ]; then
    vrm_file=/pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac.nc
  else
    echo "unknown usecase: ${usecase}"
    exit 1
  fi

  ${cdo} -selvar,veg_ratio_max ${vrm_file} vrm.nc
  ${cdo} -selvar,land ${vrm_file} land.nc
  ${cdo} -selvar,fract_veg ${vrm_file} fract_veg.nc
  ${cdo} -mul land.nc fract_veg.nc veg.nc
  ${cdo} -mul veg.nc vrm.nc veg_vrm.nc
  rm land.nc fract_veg.nc veg.nc vrm.nc

  # divide tmp_remapped.nc by veg_vrm.nc for all variables:
  # - when veg_vrm is zero, then set the result to zero instead of NaN
  # - retain NaN when the original value is NaN
  # - when the result is >1, set it to 1 (this concerns <300 grid cells with very low veg_vrm in JSB and
  #   relatively high transition rates in the LUH3 data)
  varList=$(${cdo} -showvar tmp_remapped.nc)
  years=$(seq ${start_year} ${end_year})
  for year in ${years}; do
    for var in ${varList}; do
      ${cdo} -selvar,${var} -selyear,${year} tmp_remapped.nc tmp_year.nc
      ${cdo} -ifthenelse veg_vrm.nc -div tmp_year.nc veg_vrm.nc veg_vrm.nc tmp2_year.nc
      ${cdo} -expr,"${var}=(${var}>1)?1:${var}" tmp2_year.nc tmp_scaled_${year}_${var}.nc
      rm tmp_year.nc tmp2_year.nc
    done
    # JSBACH and the JSBACH scripting expect one file per year
    if [ -f LUH_transitions_r2b4_${usecase}_${year}.nc ]; then
      rm LUH_transitions_r2b4_${usecase}_${year}.nc
    fi
    ${cdo} -merge tmp_scaled_${year}_*.nc LUH_transitions_r2b4_${usecase}_${year}.nc
    rm tmp_scaled_${year}_*.nc
  done
  rm tmp_remapped.nc veg_vrm.nc

  # reduce time dimension
  files=$( ls LUH_transitions_r2b4_${usecase}_????.nc )

  for file in ${files}; do
    ${cdo} --reduce_dim copy ${file} tmp.nc
    mv tmp.nc ${file}
  done

else
  echo "No remapping to ICON grid..."
  echo "No scaling with vrm..."

  # JSBACH and the JSBACH scripting expect one file per year
  ${cdo} -splityear LUH_trans_${start_year}-${end_year}.nc LUH_transitions_${usecase}_
  rm LUH_trans_${start_year}-${end_year}.nc

  # reduce time dimension
  files=$( ls LUH_transitions_${usecase}_????.nc )

  for file in ${files}; do
    ${cdo} --reduce_dim copy ${file} tmp.nc
    mv tmp.nc ${file}
  done
fi