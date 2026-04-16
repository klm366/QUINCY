#!/usr/bin/env bash

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

#----------------------------------------------------------------------
# Collects those QUINCY vegetation, hydrology, and soil snow energy variables from the output nc files
# to generate the soil forcing files for the  QUINCY SOIL model
#
#
# Notes: 
# - The collected variables should be in timestep of PT30M of one year period
# - They should be in sync with the variables read in update_read_soil_forcing
#
# Inspired by cpool_file_from_restart_files.sh

module unload cdo
module unload nco
module load nco/5.0.6-gcc-11.2.0
module load cdo/2.5.0-gcc-11.2.0

cdo="cdo -s -P 8"
#----------------------------------------------------------------------
# paths and filenames

input_path="./"
input_file_sp="qsoil_forcing_sp.nc"
input_file1_litter="qsoil_forcing_litter1.nc"
input_file2_litter="qsoil_forcing_litter2.nc"
output_path="./"
output_file_sp="bc_qsoil_physics.nc"
output_file_litter="bc_qsoil_litterfall.nc"

# usecase dependent number of pfts
nr_pft=13

# the name of the extra dimension in the icon output
name_of_soildim=layers_5

# selected soil model (simple_1d or jsm)
soilModel=simple_1d

# Assert that the output file does not already exist
if [ -f  ${output_path}/${output_file_sp} ]; then
  echo "File ${output_path}/${output_file_sp} already exists. Please delete or rename to proceed"
  rm -f ${output_path}/${output_file_sp}
  # exit
fi

if [ -f  ${output_path}/${output_file_litter} ]; then
  echo "File ${output_path}/${output_file_litter} already exists. Please delete or rename to proceed"
  rm -f ${output_path}/${output_file_litter}
  # exit
fi
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# To derive the variable names as used in the output file
# - for the litterfall bgcms
#   - all required elements
#   - all required compartments 
#   - bgcm main prefix (veg and seedbed litterfall, respectively)
# - for 'normal' variables
#   - process prefix (currently veg, sse or hydro)
#   - variable names
# - for all variables
#   - tile name postfix (currently only pfts)

#------------------
#-- elements
elements="C N P C13 C14 N15"

#------------------
#-- vegetation, veg litter, seedbed litter
# bgcm
veg_bgcm_prefix="veg_veg_bgcm_pool"
veg_litter_bgcm_prefix="veg_veg_bgcm_litterfall"
seedbed_litter_bgcm_prefix="veg_veg_bgcm_seed_bed_litterfall"
veg_compartments_1="leaf  fine_root  coarse_root  sap_wood"
veg_compartments_2="heart_wood  labile  reserve  fruit"

#------------------
# other forcing variables
veg_variable_list="veg_root_fraction_sl veg_demand_uptake_n_tuptake_mavg veg_demand_uptake_p_tuptake_mavg"
variable_list="hydro_frac_wtr_transp_down_sl  hydro_frac_w_lat_loss_sl
                   hydro_wtr_soil_sl  hydro_vol_porosity_sl  sse_t_soil_sl  a2l_nhx_deposition
                   a2l_noy_deposition  a2l_nhx_n15_deposition a2l_noy_n15_deposition
                   a2l_p_deposition"

#------------------
#-- tiles: all pfts
tile_list=""
i=1
while [ ${i} -le ${nr_pft} ]; do
  if [ ${i} -lt 10 ]; then
    tile_list="${tile_list} pft0${i}"
  else
    tile_list="${tile_list} pft${i}"
  fi
  (( i = i + 1 ))
done

#----------------------------------------------------------------------
which ncatted || {
. ${MODULESHOME}/init/bash
  module load nco
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Generate list of to be extracted variables for output_file_litter
var_list_sp=""
var_list_litter1=""
var_list_litter2=""

# veg fine root C
for tile in ${tile_list} ; do
  var_list_sp="${var_list_sp},${veg_bgcm_prefix}_fine_root_C_${tile}"
done

# veg bgcm litterfall
for element in ${elements}; do
  for tile in ${tile_list} ; do
    for compartment in ${veg_compartments_1} ; do
      var_list_litter1="${var_list_litter1},${veg_litter_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  
    for compartment in ${veg_compartments_2} ; do
      var_list_litter2="${var_list_litter2},${veg_litter_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  done
done

# seedbed bgcm litterfall
for element in ${elements}; do
  for tile in ${tile_list} ; do
    var_list_litter2="${var_list_litter2},${seedbed_litter_bgcm_prefix}_${element}_${tile}"
  done
done

# extract the variables from restart files
tmp_litter_1="${output_path}/tmp_litter_1.nc"


${cdo} --no_history --reduce_dim -selvar${var_list_sp} ${input_path}/${input_file_sp} ${output_path}/${tmp_litter_1}

# # merge tmp files with cdo
# ${cdo} --no_history -b 64 merge "$tmp_litter_1" "${input_path}/${input_file1_litter}" "${input_path}/${input_file2_litter}" "${output_path}/${output_file_litter}"

# merge tmp files with nco, optional
cp -f "${output_path}/${tmp_litter_1}" "${output_path}/${output_file_litter}"
ncks -A "${input_path}/${input_file1_litter}" "${output_path}/${output_file_litter}"
ncks -A "${input_path}/${input_file2_litter}" "${output_path}/${output_file_litter}"

# rename dimension
ncrename -h -d .ncells,cell ${output_path}/${output_file_litter}
# remove tmp files
rm "$tmp_litter_1"

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Generate list of to be extracted variables for output_file_sp
var_list_water=""
var_list_energy=""

# veg variables
for var in ${veg_variable_list} ; do
  for tile in ${tile_list} ; do
    var_list_water="${var_list_water},${var}_${tile}"
  done
done

# other variables
for var in ${variable_list}; do
  v="${var}_box"

  # Inspect dimensions of v in the INPUT file
  # (works well with NCO; adjust if your environment differs)
  dims=$(ncks -m -v "$v" "${input_path}/${input_file_sp}" 2>/dev/null | head -n 40)

  if echo "$dims" | grep -q "soil_depth_water"; then
    var_list_water="${var_list_water},${v}"
  else
    var_list_energy="${var_list_energy},${v}"
  fi
done

# extract the variables from restart file
tmp_water="${output_path}/tmp_sp_water.nc"
tmp_energy="${output_path}/tmp_sp_energy.nc"

# helper: extract only if list has something after the leading comma
extract_if_nonempty () {
  local list="$1" infile="$2" outfile="$3"
  if [[ -n "${list#,}" ]]; then
    ${cdo} --no_history -b 64 --reduce_dim -selvar${list} "$infile" "$outfile"
  fi
}

extract_if_nonempty "$var_list_water"  "${input_path}/${input_file_sp}" "$tmp_water"
extract_if_nonempty "$var_list_energy" "${input_path}/${input_file_sp}" "$tmp_energy"

# Horizontal dim
for f in "$tmp_water" "$tmp_energy"; do
  [[ -f "$f" ]] && ncrename -h -d .ncells,cell "$f"
done

# Soil dim -> soillev (only inside the matching file)
[[ -f "$tmp_water"  ]] && ncrename -h -d .soil_depth_water,soillev  "$tmp_water"
[[ -f "$tmp_energy" ]] && ncrename -h -d .soil_depth_energy,soillev "$tmp_energy"

${cdo} --no_history -b 64 merge "$tmp_water" "$tmp_energy" "${output_path}/${output_file_sp}"
# remove tmp files
rm "$tmp_water" "$tmp_energy"


#------------------------------------------------------------------------------
# set global attributes
ncatted -h -O -a ,global,d,, ${output_path}/${output_file_sp}
ncatted -h -O -a comment,global,o,c,"Soil physics forcing extracted from output file ${input_file}" ${output_path}/${output_file_sp}

ncatted -h -O -a ,global,d,, ${output_path}/${output_file_litter}
ncatted -h -O -a comment,global,o,c,"Litterfall forcing extracted from output file ${input_file}" ${output_path}/${output_file_litter}