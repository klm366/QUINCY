#!/bin/ksh

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
set -e
#-----------------------------------------------------------------------------
#  Generate
#     ICON-Land land cover transition files
#
# The script is part of a script series to generate ICON-Land initial files
# started by master script "create_jsbach_ini_files.sh".
#-----------------------------------------------------------------------------

# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}       # Extpar file name

icon_grid=${icon_grid}
atmGridID=${atmGridID}
remap_scheme=${remap_scheme}

path_bc=${bc_file_dir}
year_list=${year_list}
work_dir=${work_dir}
pool_prepare=${pool_prepare}

cdo=${cdo}
rm=$rm
cp=$cp

# Local variables
prog=$(basename $0)
pft_tag="12pfts"  # recommended: 12pfts
                  # Note that "11pfts" is also possible for transitions, but consistent cover fractions
                  # (i.e. based on LUH3 states) are not available for 11 PFTs at the moment.
file_bc_frac=${path_bc}/bc_land_frac.nc # generated in generate_bc_land_fractions.sh
file_bc_trans=bc_land_trans
clean_up=true  # Remove intermediate files

# Input files for transitions:
#       aggregated transition rates based on LUH3 data generated with luh3/aggregate_LUH_transitions.sh
#luh3_trans_root=${pool_prepare}/025/LUH3/r0001/LUH_transitions   # without _${pft_tag}_${year}.nc
luh3_trans_root=/work/bm0891/m300218/LUH3/prepare_transitions/LUH_transitions # without _${pft_tag}_${year}.nc


start_year=$(echo ${year_list} | cut -f1 -d" ")
year=${start_year}

prep_025_wgts_file=${luh3_trans_root}_${pft_tag}_${year}.nc

cd ${work_dir}

echo "${prog}: Generating land cover transition files for JSBACH setup."

if [[ ! -f ${file_bc_frac} ]]; then
  echo "ERROR: ${prog}: bc fraction file ${file_bc_frac} not found."
  echo "       Please run generate_bc_land_fractions.sh first."
  exit 1
fi

if [[ ! -f non-glac-land-mask.nc ]]; then
  echo "ERROR: ${prog}: non-glac-land-mask.nc not found."
  echo "       Please run bc_files_from_extpar.sh first."
  exit 1
fi

# -----------------------------------
#  Land cover transitions
# -----------------------------------
# ------------------
# remapping to ICON grid and rescaling
# - Remapping from 0.25 degree lat-lon grid to ICON grid
# - In LUH3, land-use transitions are given as the fraction of each grid-cell that transitions from one land-use
#   state to another in a given year. Bare lands (e.g. deserts) are included in natural vegetation in LUH3, but
#   in JSBACH, bare soil is implicit. Therefore, the transition rates from LUH3 need to be scaled with veg_ratio_max.
# ------------------

# -----------------------------------
#  1. find variables for scaling
# -----------------------------------
# veg_ratio_max is constant in time in the bc files; we only need to extract it once.
[[ -f vrm.nc ]] && ${rm} vrm.nc
${cdo} -selvar,veg_ratio_max ${file_bc_frac} vrm.nc
var_list=$(${cdo} -showvar ${luh3_trans_root}_${pft_tag}_${year}.nc)

if [[ ! -f wgt025.nc ]]; then
  ${cdo} gen${remap_scheme},${icon_grid} -setmisstoc,0. ${prep_025_wgts_file} wgt025.nc
fi

for year in ${year_list}; do

  if [[ ! -f ${luh3_trans_root}_${pft_tag}_${year}.nc ]]; then
    echo "ERROR: ${prog}: file with aggregated transition rates ${luh3_trans_root}_${pft_tag}_${year}.nc not found."
    echo "       Please run luh3/aggregate_LUH_transitions.sh first."
    exit 1
  fi

  # -----------------------------------
  #  2. remapping
  # -----------------------------------

  ${cdo} mul -remap,${icon_grid},wgt025.nc -setmisstodis ${luh3_trans_root}_${pft_tag}_${year}.nc \
         non-glac-land-mask.nc ${atmGridID}_trans_${year}.nc

  # -----------------------------------
  #  3. scaling of transitions
  # -----------------------------------

  # In JSBACH bare soil is implicit. Therefore, cover fractions are scaled up by veg_ratio_max (elsewhere)
  # and this scaling is applied (here) to the transitions as well.
  # Divide ${atmGridID}_trans_${year}.nc by veg_ratio_max for all variables:
  # - when veg_vrm is zero, then set the result to zero instead of NaN
  # - retain NaN when the original value is NaN
  # - when the result is >1, set it to 1 (this applies to a number of grid cells with very low veg_ratio_max in JSB and
  #   relatively high transition rates in the LUH3 dataset)

  for var in ${var_list}; do
    {
    ${cdo} -selvar,${var} ${atmGridID}_trans_${year}.nc tmp_${var}.nc
    ${cdo} -ifthenelse vrm.nc -div tmp_${var}.nc vrm.nc vrm.nc tmp2_${var}.nc
    ${cdo} -expr,"${var}=(${var}>1)?1:${var}" tmp2_${var}.nc tmp_scaled_${year}_${var}.nc
    ${rm} tmp_${var}.nc tmp2_${var}.nc
    } &  # comment out '&' for serial processing
  done
  wait

  # JSBACH and the JSBACH scripting expect one file per year
  if [[ -f LUH_transitions_${atmGridID}_${pft_tag}_${year}.nc ]]; then
    ${rm} LUH_transitions_${atmGridID}_${pft_tag}_${year}.nc
  fi
  ${cdo} -merge tmp_scaled_${year}_*.nc tmp.nc
  ${cdo} --reduce_dim copy tmp.nc ${file_bc_trans}_${pft_tag}_${year}.nc
  ${cdo} --no_history setattribute,history="${history_att}" ${file_bc_trans}_${pft_tag}_${year}.nc \
       ${path_bc}/${file_bc_trans}_${pft_tag}_${year}.nc
  ncks -A -x -h ${luh3_trans_root}_${pft_tag}_${year}.nc ${path_bc}/${file_bc_trans}_${pft_tag}_${year}.nc
  ${rm} tmp_scaled_${year}_*.nc ${atmGridID}_trans_${year}.nc
  ${rm} tmp.nc ${file_bc_trans}_${pft_tag}_${year}.nc

  echo "${prog}:     ${file_bc_trans}_${pft_tag}_${year}.nc             done"
done

${rm} vrm.nc

echo "${prog}:  done"

exit 0