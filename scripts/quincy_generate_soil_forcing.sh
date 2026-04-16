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
# Collect those QUINCY vegetation, hydrology, and soil snow energy variables
# from ICON-Land output files in order to generate forcing files for the
# QUINCY SOIL model.
#
#
# Notes:
# - This script supports two new source layouts:
#   1) one yearly output file with daily records
#   2) 365 daily output files with PT30M records
# - The collected variables must stay in sync with the variables read in
#   update_read_soil_forcing.
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
output_path="./"
output_file_sp="bc_qsoil_physics.nc"
output_file_litter="bc_qsoil_litterfall.nc"
output_file_sp_prefix="bc_qsoil_physics"
output_file_litter_prefix="bc_qsoil_litterfall"

# Source layout:
# - yearly_p1d : one yearly file with daily records
# - daily_pt30m: 365 daily files with PT30M records
source_layout="yearly_p1d"

# ICON output naming convention:
#   <exp_name>_<input_selector>_<timestamp>.nc
# e.g.
#   quincy_R2B5_land-CN-gswp3_S0_lnd_basic_ml_18500101T000000Z.nc
exp_name="quincy_R2B5_land-CN-gswp3_S0"
input_selector="lnd_basic_ml"
forcing_year=1850

# usecase dependent number of pfts
nr_pft=13

# the name of the extra dimension in the icon output
name_of_soildim=layers_5

# selected soil model (simple_1d or jsm)
soilModel=simple_1d

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
# Generate lists of to be extracted variables
var_list_sp=""
var_list_litter=""

# veg fine root C
for tile in ${tile_list} ; do
  var_list_litter="${var_list_litter},${veg_bgcm_prefix}_fine_root_C_${tile}"
done

# veg bgcm litterfall
for element in ${elements}; do
  for tile in ${tile_list} ; do
    for compartment in ${veg_compartments_1} ; do
      var_list_litter="${var_list_litter},${veg_litter_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  
    for compartment in ${veg_compartments_2} ; do
      var_list_litter="${var_list_litter},${veg_litter_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  done
done

# seedbed bgcm litterfall
for element in ${elements}; do
  for tile in ${tile_list} ; do
    var_list_litter="${var_list_litter},${seedbed_litter_bgcm_prefix}_${element}_${tile}"
  done
done
var_list_water=""
var_list_energy=""

# veg variables
for var in ${veg_variable_list} ; do
  for tile in ${tile_list} ; do
    var_list_water="${var_list_water},${var}_${tile}"
  done
done

# other variables
extract_if_nonempty () {
  local list="$1" infile="$2" outfile="$3"
  if [[ -n "${list#,}" ]]; then
    ${cdo} --no_history -b 64 --reduce_dim -selvar${list} "$infile" "$outfile"
  fi
}

timestamp_for_day () {
  local year="$1"
  local day_of_year="$2"
  date -u -d "${year}-01-01 +$((day_of_year - 1)) day" +%Y%m%dT%H%M%SZ
}

doy_for_day () {
  local day_of_year="$1"
  printf '%03d' "${day_of_year}"
}

build_source_file_list () {
  local year="$1"
  local source_files=()
  local timestamp
  local input_file
  local day

  case "${source_layout}" in
    yearly_p1d)
      timestamp=$(timestamp_for_day "${year}" 1)
      input_file="${input_path}/${exp_name}_${input_selector}_${timestamp}.nc"
      source_files+=("${input_file}")
      ;;
    daily_pt30m)
      for day in $(seq 1 365); do
        timestamp=$(timestamp_for_day "${year}" "${day}")
        input_file="${input_path}/${exp_name}_${input_selector}_${timestamp}.nc"
        source_files+=("${input_file}")
      done
      ;;
    *)
      echo "Unsupported source_layout: ${source_layout}"
      exit 1
      ;;
  esac

  printf '%s\n' "${source_files[@]}"
}

prepare_var_groups () {
  local input_file="$1"
  local dims
  local var
  local v

  var_list_water=""
  var_list_energy=""

  for var in ${veg_variable_list} ; do
    for tile in ${tile_list} ; do
      var_list_water="${var_list_water},${var}_${tile}"
    done
  done

  for var in ${variable_list}; do
    v="${var}_box"
    dims=$(ncks -m -v "$v" "${input_file}" 2>/dev/null | head -n 40)

    if echo "$dims" | grep -q "soil_depth_water"; then
      var_list_water="${var_list_water},${v}"
    else
      var_list_energy="${var_list_energy},${v}"
    fi
  done
}

merge_sp_files () {
  local water_file="$1"
  local energy_file="$2"
  local output_file="$3"

  if [[ -f "${water_file}" && -f "${energy_file}" ]]; then
    ${cdo} --no_history -b 64 merge "${water_file}" "${energy_file}" "${output_file}"
  elif [[ -f "${water_file}" ]]; then
    cp -f "${water_file}" "${output_file}"
  elif [[ -f "${energy_file}" ]]; then
    cp -f "${energy_file}" "${output_file}"
  else
    echo "No soil-physics variables were extracted for ${output_file}"
    exit 1
  fi
}

add_global_comment () {
  local output_file="$1"
  local comment="$2"

  ncatted -h -O -a ,global,d,, "${output_file}"
  ncatted -h -O -a comment,global,o,c,"${comment}" "${output_file}"
}

process_source_file () {
  local input_file="$1"
  local output_sp_file="$2"
  local output_litter_file="$3"
  local tmp_water="${output_path}/tmp_sp_water.nc"
  local tmp_energy="${output_path}/tmp_sp_energy.nc"
  local tmp_litter="${output_path}/tmp_litter.nc"
  local f

  if [[ ! -f "${input_file}" ]]; then
    echo "Missing input file: ${input_file}"
    exit 1
  fi

  rm -f "${output_sp_file}" "${output_litter_file}" "${tmp_water}" "${tmp_energy}" "${tmp_litter}"

  prepare_var_groups "${input_file}"

  extract_if_nonempty "${var_list_litter}" "${input_file}" "${tmp_litter}"
  extract_if_nonempty "${var_list_water}"  "${input_file}" "${tmp_water}"
  extract_if_nonempty "${var_list_energy}" "${input_file}" "${tmp_energy}"

  if [[ -f "${tmp_litter}" ]]; then
    ncrename -h -d .ncells,cell "${tmp_litter}"
    cp -f "${tmp_litter}" "${output_litter_file}"
  else
    echo "No litter variables were extracted for ${input_file}"
    exit 1
  fi

  for f in "${tmp_water}" "${tmp_energy}"; do
    [[ -f "${f}" ]] && ncrename -h -d .ncells,cell "${f}"
  done

  [[ -f "${tmp_water}"  ]] && ncrename -h -d .soil_depth_water,soillev  "${tmp_water}"
  [[ -f "${tmp_energy}" ]] && ncrename -h -d .soil_depth_energy,soillev "${tmp_energy}"

  merge_sp_files "${tmp_water}" "${tmp_energy}" "${output_sp_file}"

  add_global_comment "${output_sp_file}" "Soil physics forcing extracted from output file ${input_file}"
  add_global_comment "${output_litter_file}" "Litterfall forcing extracted from output file ${input_file}"

  rm -f "${tmp_water}" "${tmp_energy}" "${tmp_litter}"
}

main () {
  local source_files=()
  local input_file
  local day_index
  local doy_string
  local output_sp_file
  local output_litter_file

  mapfile -t source_files < <(build_source_file_list "${forcing_year}")

  case "${source_layout}" in
    yearly_p1d)
      process_source_file "${source_files[0]}" \
        "${output_path}/${output_file_sp}" \
        "${output_path}/${output_file_litter}"
      ;;
    daily_pt30m)
      day_index=0
      for input_file in "${source_files[@]}"; do
        (( day_index = day_index + 1 ))
        doy_string=$(doy_for_day "${day_index}")

        output_sp_file="${output_path}/${output_file_sp_prefix}_${doy_string}.nc"
        output_litter_file="${output_path}/${output_file_litter_prefix}_${doy_string}.nc"

        process_source_file "${input_file}" "${output_sp_file}" "${output_litter_file}"
      done
      ;;
  esac
}

main "$@"
