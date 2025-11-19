#!/bin/ksh

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------
set -e
#-----------------------------------------------------------------------------
#   Generate
#     ICON-Land soil initial conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh".
#
#   ic_land_soil:
#   - Initial soil mositure and snow cover

# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}

icon_grid=${icon_grid}
atmGridID=${atmGridID}

path_bc=${bc_file_dir}
work_dir=${work_dir}
pool_prepare=${pool_prepare}

cdo=${cdo}
rm=$rm
cp=$cp

# Local variables
prog=$(basename $0)
file_ic_soil=ic_land_soil
clean_up=true  # Remove intermediate files

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Soil initial conditions
# -----------------------------
echo "${prog}: Generating ${file_ic_soil}:"

# The order of variables should not be changed
var_list="init_moist surf_temp layer_moist snow"

# Preparations
jan_surf_file=${pool_prepare}/05/ECHAM6/05_jan_surf.nc
ini_mon=01
${cdo}  ifthen land.nc -mul notglac.nc -remap,${icon_grid},wgt05.nc -setmisstodis \
  -ifthen -selvar,SLM ${jan_surf_file} -selvar,WS ${jan_surf_file}  ${atmGridID}_WS.tmp

# init_moist
# ----------------------
# Initial soil moisture is limited by the field capacity and should be higher than
# 50% to prevent vegetation die back due to initialization.
# In desert regions (veg_ratio_max < 20%), soil moisture is initialized below the
# wilting point.
${cdo} mul ${atmGridID}_soil_field_cap.keep ${atmGridID}_root_depth.keep \
   field_cap_rootzone.tmp
${cdo} ltc,0.2 ${atmGridID}_veg_ratio_max_12pfts.keep desert-mask.nc
${cdo} ifthenelse desert-mask.nc \
         -mulc,0.3 field_cap_rootzone.tmp \
         -mulc,0.5 field_cap_rootzone.tmp  init_moist_min.tmp
${cdo} -setattribute,init_moist@long_name="Soil water content" \
    -setattribute,init_moist@units="m water equivalent" \
    -setmissval,-9.e+33 -setname,init_moist \
    -ifthen land.nc \
    -min field_cap_rootzone.tmp -max ${atmGridID}_WS.tmp init_moist_min.tmp \
    ${atmGridID}_init_moist.nc


for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    init_moist )
  # ----------------------
      # Done above - outside parallel loop
      ;;
  # ----------------------
    layer_moist )
  # ----------------------
      levels="0.065 0.319 1.232 4.134 9.834"

      # Relative root zone moisture
      ${cdo} setmisstoc,0. -div ${atmGridID}_init_moist.nc ${atmGridID}_root_depth.keep \
          init_moist_rel.tmp

      above=0.
      for lev in ${levels}; do
         ldepth=$(echo "${lev} - ${above}" | bc)                    # Current layer depth
         ${cdo} minc,${lev} -subc,${above} ${atmGridID}_soil_depth.keep sdepth.tmp  # Soil depth in layer
         ${cdo} ifthenelse -gtc,0 sdepth.tmp \
               -mul init_moist_rel.tmp sdepth.tmp \
               zero.nc   ${var}_${lev}.tmp1
         ncecat -u soillev ${var}_${lev}.tmp1 ${var}_${lev}.tmp2    # Add new record dimension "soillev"
         ncecat ${var}_${lev}.tmp2 ${var}_${lev}.tmp3               # Convert soilev to non-record dim.
         ncwa -a record ${var}_${lev}.tmp3 ${var}_${lev}.tmp4       # Remove degenerate record dim.
         ${cdo} setlevel,${lev} ${var}_${lev}.tmp4 ${var}_${lev}.tmp  # Set level
         above=${lev}
      done
      $rm -f ${atmGridID}_${var}.nc
      ${cdo} -setattribute,${var}@long_name="Soil water content" \
          -setattribute,${var}@units="m" -setname,${var} \
          -setattribute,soillev@long_name="Soil layer (lower boundary)" \
          -setattribute,soillev@units="m" \
          -setmissval,-9.e+33 \
          -ifthen land.nc \
          -merge ${var}_?????.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    snow )
  # ----------------------
      # Limit to 20 cm (water equivalent)
      ${cdo} -setattribute,${var}@long_name="Snow depth" -setattribute,${var}@units="m water equivalent" \
          -setname,${var} -setmissval,-9.e+33 -ifthen land.nc -mul non-glac-land-mask.nc -minc,0.2 -mulc,0.001 \
          -seldate,1111-${ini_mon}-11 -selvar,W_SNOW ${extpar_file} \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    surf_temp )
  # ----------------------
      ${cdo} -setname,${var} -setmissval,-9.e+33 -selvar,T_2M_CLIM ${extpar_file} \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_soil file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait

# Merge newly calculated variables to the ic soil file
# -------------
[[ -f ${file_ic_soil}.tmp ]] && ${rm} ${file_ic_soil}.tmp
for var in ${var_list}; do
  if [[ ! -f ${file_ic_soil}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_ic_soil}.tmp
  else
    # following variable
    mv ${file_ic_soil}.tmp ${file_ic_soil}.tmp2
    ${cdo} -O merge ${file_ic_soil}.tmp2  ${atmGridID}_${var}.nc ${file_ic_soil}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_ic_soil}.tmp \
       ${path_bc}/${file_ic_soil}.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
fi
echo "${prog}:     ${file_ic_soil}.nc          done"

exit 0
