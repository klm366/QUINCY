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

#------------------------------------------------------------------------------
# Generate a static analytical steady-state SB initialisation NetCDF file from
# an ICON-QUINCY output file with a time dimension.
#
# The resulting file is meant to be used directly as:
#   sb_analytical_init_filename
#
# Workflow:
# - take time means for all rate / modifier / state variables
# - aggregate PFT-resolved litter-formation C fluxes to one grid-cell value
# - convert litter formation from [mol m-3 timestep-1] to [mol m-3 s-1]
# - rename variables to the exact names expected by sb_read_analytical_init_vars
#
# Notes:
# - The ICON output file provides the three litter-formation C variables only
#   per PFT, so they are summed across PFTs here.
# - All other analytical-init variables are assumed to be available directly
#   in the ICON output stream on the box tile.
#------------------------------------------------------------------------------

set -euo pipefail

module unload cdo >/dev/null 2>&1 || true
module unload nco >/dev/null 2>&1 || true
module load nco/5.0.6-gcc-11.2.0
module load cdo/2.5.0-gcc-11.2.0

cdo="cdo -s -P 8"

#------------------------------------------------------------------------------
# paths and filenames

input_path="./"
input_file="icon_sb_analytical_init_source.nc"
output_path="./"
output_file="bc_land_sb_analytical_init.nc"

# Number of PFTs in the ICON output stream
nr_pft=13

# Model timestep used for the litter-formation outputs written in
# [mol m-3 timestep-1]. Adjust if the source simulation used another dtime.
dtime_seconds=1800

input_filepath="${input_path%/}/${input_file}"
output_filepath="${output_path%/}/${output_file}"

if [[ ! -f "${input_filepath}" ]]; then
  echo "Input file does not exist: ${input_filepath}"
  exit 1
fi

if [[ -f "${output_filepath}" ]]; then
  echo "Output file already exists: ${output_filepath}"
  exit 1
fi

tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

tmp_mean="${tmp_dir}/analytical_init_timemean.nc"
tmp_litter="${tmp_dir}/analytical_init_litter.nc"
tmp_all="${tmp_dir}/analytical_init_all.nc"

has_var () {
  local var_name="$1"
  cdo -s showname "${input_filepath}" | tr ' ' '\n' | grep -Fxq "${var_name}"
}

make_pft_list () {
  local prefix="$1"
  local vars=()
  local i
  for (( i=1; i<=nr_pft; i++ )); do
    vars+=("${prefix}_pft$(printf '%02d' "${i}")")
  done
  printf '%s\n' "${vars[@]}"
}

required_direct_mappings=(
  "hydro_frac_wtr_transp_down_sl_box:frac_wtr_transp_down_sl"
  "hydro_frac_w_lat_loss_sl_box:frac_w_lat_loss_sl"
  "sb_rtm_decomposition_box:rtm_decomposition"
  "sb_rmm_decomposition_box:rmm_decomposition"
  "sb_rtm_depolymerisation_box:rtm_depolymerisation"
  "sb_rmm_depolymerisation_box:rmm_depolymerisation"
  "sb_rtm_sorption_box:rtm_sorption"
  "sb_rmm_sorption_box:rmm_sorption"
  "sb_rtm_desorption_box:rtm_desorption"
  "sb_rmm_desorption_box:rmm_desorption"
  "sb_rtm_mic_uptake_box:rtm_mic_uptake"
  "sb_rmm_mic_uptake_box:rmm_mic_uptake"
  "sb_rtm_hsc_box:rtm_hsc"
  "sb_rmm_hsc_box:rmm_hsc"
  "sb_enzyme_frac_AP_box:enzyme_frac_ap"
  "sse_bulk_dens_sl_box:bulk_dens_sl"
  "sse_silt_sl_box:silt_sl"
  "sse_clay_sl_box:clay_sl"
  "sse_t_soil_sl_box:t_soil_sl"
  "hydro_wtr_soil_sl_box:wtr_soil_sl"
  "hydro_wtr_soil_pot_sl_box:wtr_soil_pot_sl"
  "hydro_ice_soil_sl_box:ice_soil_sl"
  "hydro_vol_field_cap_sl_box:vol_field_cap_sl"
  "hydro_soil_depth_sl_box:soil_depth_sl"
)

litter_mappings=(
  "sb_sb_bgcm_formation_sol_litter_C:formation_soluable_litter_c"
  "sb_sb_bgcm_formation_pol_litter_C:formation_polymeric_litter_c"
  "sb_sb_bgcm_formation_woo_litter_C:formation_woody_litter_c"
)

selected_direct_vars=()
direct_renames=()

for mapping in "${required_direct_mappings[@]}"; do
  src_var="${mapping%%:*}"
  dst_var="${mapping##*:}"
  if ! has_var "${src_var}"; then
    echo "Missing required variable in ICON output: ${src_var}"
    exit 1
  fi
  selected_direct_vars+=("${src_var}")
  direct_renames+=("${src_var}:${dst_var}")
done

selected_litter_vars=()
for mapping in "${litter_mappings[@]}"; do
  src_prefix="${mapping%%:*}"
  while IFS= read -r var_name; do
    if ! has_var "${var_name}"; then
      echo "Missing required PFT litter variable in ICON output: ${var_name}"
      exit 1
    fi
    selected_litter_vars+=("${var_name}")
  done < <(make_pft_list "${src_prefix}")
done

if [[ ${#selected_direct_vars[@]} -gt 0 ]]; then
  selvar_direct="$(IFS=,; echo "${selected_direct_vars[*]}")"
  ${cdo} -b 64 timmean -selvar,"${selvar_direct}" "${input_filepath}" "${tmp_mean}"
else
  ncks -O -x -v time "${input_filepath}" "${tmp_mean}"
fi

selvar_litter="$(IFS=,; echo "${selected_litter_vars[*]}")"
${cdo} -b 64 timmean -selvar,"${selvar_litter}" "${input_filepath}" "${tmp_litter}"

for mapping in "${litter_mappings[@]}"; do
  src_prefix="${mapping%%:*}"
  dst_var="${mapping##*:}"

  expr="${dst_var}=("
  first_term=1
  while IFS= read -r var_name; do
    if [[ ${first_term} -eq 1 ]]; then
      expr+="${var_name}"
      first_term=0
    else
      expr+="+${var_name}"
    fi
  done < <(make_pft_list "${src_prefix}")
  expr+=")/${dtime_seconds}.0"

  ncap2 -O -s "${expr}" "${tmp_litter}" "${tmp_litter}"
done

# Keep only the aggregated litter fields.
tmp_litter_keep=()
for mapping in "${litter_mappings[@]}"; do
  tmp_litter_keep+=("${mapping##*:}")
done
selvar_litter_keep="$(IFS=,; echo "${tmp_litter_keep[*]}")"
ncks -O -v "${selvar_litter_keep}" "${tmp_litter}" "${tmp_litter}"

# Remove singleton time dimension after timmean.
if ncks -m "${tmp_mean}" 2>/dev/null | grep -Eq 'time[[:space:]]*=|UNLIMITED[[:space:]]*//[[:space:]]*\(time'; then
  ncwa -O -a time "${tmp_mean}" "${tmp_mean}"
fi
if ncks -m "${tmp_litter}" 2>/dev/null | grep -Eq 'time[[:space:]]*=|UNLIMITED[[:space:]]*//[[:space:]]*\(time'; then
  ncwa -O -a time "${tmp_litter}" "${tmp_litter}"
fi

for mapping in "${direct_renames[@]}"; do
  src_var="${mapping%%:*}"
  dst_var="${mapping##*:}"
  if [[ "${src_var}" != "${dst_var}" ]]; then
    ncrename -O -v "${src_var},${dst_var}" "${tmp_mean}"
  fi
done

# Normalise dimensional names as far as possible for the jsbach reader.
if ncks -m "${tmp_mean}" 2>/dev/null | grep -Fq "soil_depth_water"; then
  ncrename -O -d "soil_depth_water,soillev" "${tmp_mean}"
fi
if ncks -m "${tmp_litter}" 2>/dev/null | grep -Fq "soil_depth_water"; then
  ncrename -O -d "soil_depth_water,soillev" "${tmp_litter}"
fi

# Merge the direct analytical variables and the aggregated litter inputs.
cp "${tmp_mean}" "${tmp_all}"
ncks -A "${tmp_litter}" "${tmp_all}"

ncatted -h -O -a ,global,d,, "${tmp_all}"
ncatted -h -O -a comment,global,o,c,"Time-aggregated analytical steady-state SB initialisation generated from ${input_filepath}" "${tmp_all}"
ncatted -h -O -a dtime_seconds_used,global,o,d,"${dtime_seconds}" "${tmp_all}"

mv "${tmp_all}" "${output_filepath}"

echo "Wrote ${output_filepath}"
