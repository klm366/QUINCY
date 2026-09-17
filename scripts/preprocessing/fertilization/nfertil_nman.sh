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

#!/usr/bin/env bash
#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=nfertil_nman.o%j
#SBATCH --error=nfertil_nman.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive
#SBATCH --ntasks=1
set -euo pipefail

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0

# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

#----------------------------------------------------------------------

# path and name of the target output folder
FINAL_DIR=./YOURFINALOUTPUTFOLDER
OUT_DIR="n_fert_applic"
TMP="${OUT_DIR}"/tmp
REMOVE_TMP="yes" # yes/no

MANURE_DIR="/work/mj0143/icon_land_data_pool/preprocessing/manure/"
FERT_DIR="/work/mj0143//icon_land_data_pool/trendy/2026/raw-data-20260702/nfert/"
FRAC_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"
STATES="/work/mj0143/icon_land_data_pool/trendy/2026/raw-data-20260702/luh2/states_LUH2_GCB2026.nc"
LUH_LAND_USE_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/luh/LUH_land_use"

# grid
GRID_ICON="/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc"
#GRID_ICON=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc

history="Calculate fertilisation application from NMIP"
contact="szaehle@bgc-jena.mpg.de"
comment="Calculate NO3 and NH4 fertilization application on crop and pasture"


mkdir -p "${OUT_DIR}"
mkdir -p "${TMP}"
mkdir -p "${FINAL_DIR}"

# prepare LUH file
# creating frac_crop and frac_pasture file !!! czeck with Chengs file later!!!!!!
${cdo} -O -selyear,1850/2024 ${STATES} ${TMP}/states_1850_2024.nc

# calculate crop and pasture from states
${cdo} -O -expr,'crop=c3ann+c4ann+c3per+c4per+c3nfx' \
  "${TMP}/states_1850_2024.nc" "${TMP}/crop.nc"
${cdo} -O -selname,pastr "${TMP}/states_1850_2024.nc" "${TMP}/pasture.nc"

${cdo} -O setmisstoc,0 "${TMP}/crop.nc" "${TMP}/crop0.nc"
${cdo} -O setmisstoc,0 "${TMP}/pasture.nc" "${TMP}/pasture0.nc"

${cdo} -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/crop0.nc" "${TMP}/frac_crop.nc"
${cdo} -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/pasture0.nc" "${TMP}/frac_pas.nc"

ncrename -O -h -v crop,frac_crop "${TMP}/frac_crop.nc"
ncrename -O -h -v pastr,frac_pas "${TMP}/frac_pas.nc"

${ncatted} -O -a units,frac_crop,o,c,"1" "${TMP}/frac_crop.nc"
${ncatted} -O -a units,frac_pas,o,c,"1" "${TMP}/frac_pas.nc"

###########################################################################

# scale manure by crop and pasture fractions (frac_crop, frac_pas)
${cdo} -mul "${MANURE_DIR}/Nmanure_Crop_App.nc" "${TMP}/frac_crop.nc" "${TMP}/manure_crops_raw.nc"
${cdo} -mul "${MANURE_DIR}/Nmanure_Pas_App.nc" "${TMP}/frac_pas.nc" "${TMP}/manure_pasture_raw.nc"

# get fertilization on the crop and pasture
${cdo} -mul "${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NH4.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Crop_NO3.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NO3.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Pas_NH4.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NH4.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Pas_NO3.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NO3.nc"

# conversion from g N m-2 → kg N ha-1: *1e-3*1e4 = *10
# 
${cdo} -O -setmissval,nan "${TMP}/Nfer_Crop_NH4.nc" "${TMP}/Nfer_Crop_NH4_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/Nfer_Crop_NH4_mvn.nc" "${TMP}/Nfer_Crop_NH4_clean.nc"
${cdo} -O -mulc,10 "${TMP}/Nfer_Crop_NH4_clean.nc" "${TMP}/Nfer_Crop_NH4_kgNha.nc"

${cdo} -O -setmissval,nan "${TMP}/Nfer_Crop_NO3.nc" "${TMP}/Nfer_Crop_NO3_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/Nfer_Crop_NO3_mvn.nc" "${TMP}/Nfer_Crop_NO3_clean.nc"
${cdo} -O -mulc,10 "${TMP}/Nfer_Crop_NO3_clean.nc" "${TMP}/Nfer_Crop_NO3_kgNha.nc"

${cdo} -O -setmissval,nan "${TMP}/Nfer_Pas_NH4.nc" "${TMP}/Nfer_Pas_NH4_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/Nfer_Pas_NH4_mvn.nc" "${TMP}/Nfer_Pas_NH4_clean.nc"
${cdo} -O -mulc,10 "${TMP}/Nfer_Pas_NH4_clean.nc" "${TMP}/Nfer_Pas_NH4_kgNha.nc"

${cdo} -O -setmissval,nan "${TMP}/Nfer_Pas_NO3.nc" "${TMP}/Nfer_Pas_NO3_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/Nfer_Pas_NO3_mvn.nc" "${TMP}/Nfer_Pas_NO3_clean.nc"
${cdo} -O -mulc,10 "${TMP}/Nfer_Pas_NO3_clean.nc" "${TMP}/Nfer_Pas_NO3_kgNha.nc"

${cdo} -O -setmissval,nan "${TMP}/manure_crops_raw.nc" "${TMP}/manure_crops_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/manure_crops_mvn.nc" "${TMP}/manure_crops_clean.nc"
${cdo} -O -mulc,10 "${TMP}/manure_crops_clean.nc" "${TMP}/manure_crops_kgNha.nc"

${cdo} -O -setmissval,nan "${TMP}/manure_pasture_raw.nc" "${TMP}/manure_pasture_mvn.nc"
${cdo} -O -setmisstoc,0 "${TMP}/manure_pasture_mvn.nc" "${TMP}/manure_pasture_clean.nc"
${cdo} -O -mulc,10 "${TMP}/manure_pasture_clean.nc" "${TMP}/manure_pasture_kgNha.nc"

# remove temporary files
rm -f "${TMP}/Nfer_Crop_NH4_mvn.nc" "${TMP}/Nfer_Crop_NH4_clean.nc"
rm -f "${TMP}/Nfer_Crop_NO3_mvn.nc" "${TMP}/Nfer_Crop_NO3_clean.nc"
rm -f "${TMP}/Nfer_Pas_NH4_mvn.nc" "${TMP}/Nfer_Pas_NH4_clean.nc"
rm -f "${TMP}/Nfer_Pas_NO3_mvn.nc" "${TMP}/Nfer_Pas_NO3_clean.nc"
rm -f "${TMP}/manure_crops_raw.nc" "${TMP}/manure_crops_mvn.nc" "${TMP}/manure_crops_clean.nc"
rm -f "${TMP}/manure_pasture_raw.nc" "${TMP}/manure_pasture_mvn.nc" "${TMP}/manure_pasture_clean.nc"

# remap to R2B5
${cdo} -remapcon,${GRID_ICON} "${TMP}/Nfer_Crop_NH4_kgNha.nc" "${TMP}/Nfer_Crop_NH4_kgNha_r2b5.nc"
${cdo} -remapcon,${GRID_ICON} "${TMP}/Nfer_Crop_NO3_kgNha.nc" "${TMP}/Nfer_Crop_NO3_kgNha_r2b5.nc"
${cdo} -remapcon,${GRID_ICON} "${TMP}/Nfer_Pas_NH4_kgNha.nc" "${TMP}/Nfer_Pas_NH4_kgNha_r2b5.nc"
${cdo} -remapcon,${GRID_ICON} "${TMP}/Nfer_Pas_NO3_kgNha.nc" "${TMP}/Nfer_Pas_NO3_kgNha_r2b5.nc"
${cdo} -remapcon,${GRID_ICON} "${TMP}/manure_crops_kgNha.nc" "${TMP}/manure_crops_kgNha_r2b5.nc"
${cdo} -remapcon,${GRID_ICON} "${TMP}/manure_pasture_kgNha.nc" "${TMP}/manure_pasture_kgNha_r2b5.nc"

# split into years

${cdo} splityear "${TMP}/Nfer_Crop_NH4_kgNha_r2b5.nc" "${TMP}/Nfer_Crop_NH4_"
${cdo} splityear "${TMP}/Nfer_Crop_NO3_kgNha_r2b5.nc" "${TMP}/Nfer_Crop_NO3_"
${cdo} splityear "${TMP}/Nfer_Pas_NH4_kgNha_r2b5.nc" "${TMP}/Nfer_Pas_NH4_"
${cdo} splityear "${TMP}/Nfer_Pas_NO3_kgNha_r2b5.nc" "${TMP}/Nfer_Pas_NO3_"
${cdo} splityear "${TMP}/manure_crops_kgNha_r2b5.nc" "${TMP}/manure_crops_"
${cdo} splityear "${TMP}/manure_pasture_kgNha_r2b5.nc" "${TMP}/manure_pasture_"

# going back to fert grid
for year in $(seq 1850 2024); do
  FRAC_FILE="${FRAC_DIR}/bc_land_frac_iq_13pfts_${year}.nc"
  FRAC_CROP="${TMP}/frac_crop_${year}.nc"
  FRAC_PAS="${TMP}/frac_pas_${year}.nc"

  # frac_crop = (pft11 + pft12) * fract_veg * land
  ${cdo} expr,'frac_crop=(fract_pft11+fract_pft12)*fract_veg*fract_land*notsea' \
      "${FRAC_FILE}" "${FRAC_CROP}"

  # frac_pas = (pft09 + pft10) * fract_veg * land
  ${cdo} expr,'frac_pas=(fract_pft09+fract_pft10)*fract_veg*fract_land*notsea' \
      "${FRAC_FILE}" "${FRAC_PAS}"

  # safety: frac==0 -> -1, vary small numbers >0 -> 0.001
  ${cdo} -O expr,'frac_crop=(frac_crop==0?-1:frac_crop)' \
    "${FRAC_CROP}" "${FRAC_CROP}.tmp"

  ${cdo} -O expr,'frac_crop=(frac_crop>0 && frac_crop<=0.001?0.001:frac_crop)' \
    "${FRAC_CROP}.tmp" "${FRAC_CROP}"

  ${cdo} -O expr,'frac_pas=(frac_pas==0?-1:frac_pas)' \
    "${FRAC_PAS}" "${FRAC_PAS}.tmp"

  ${cdo} -O expr,'frac_pas=(frac_pas>0 && frac_pas<=0.001?0.001:frac_pas)' \
  "${FRAC_PAS}.tmp" "${FRAC_PAS}"

  ${cdo} -O merge "${TMP}/Nfer_Crop_NH4_${year}.nc" "${FRAC_CROP}" "${TMP}/tmp1.nc"
  ${cdo} -O expr,'fertl_crops_nh4=Nfer_Crop_NH4/frac_crop' "${TMP}/tmp1.nc" "${OUT_DIR}/fertl_crops_nh4_${year}.nc"

  ${cdo} -O merge "${TMP}/Nfer_Crop_NO3_${year}.nc" "${FRAC_CROP}" "${TMP}/tmp2.nc"
  ${cdo} -O expr,'fertl_crops_no3=Nfer_Crop_NO3/frac_crop' "${TMP}/tmp2.nc" "${OUT_DIR}/fertl_crops_no3_${year}.nc"

  ${cdo} -O merge "${TMP}/Nfer_Pas_NH4_${year}.nc" "${FRAC_PAS}" "${TMP}/tmp3.nc"
  ${cdo} -O expr,'fertl_pasture_nh4=Nfer_Pas_NH4/frac_pas' "${TMP}/tmp3.nc" "${OUT_DIR}/fertl_pasture_nh4_${year}.nc"

  ${cdo} -O merge "${TMP}/Nfer_Pas_NO3_${year}.nc" "${FRAC_PAS}" "${TMP}/tmp4.nc"
  ${cdo} -O expr,'fertl_pasture_no3=Nfer_Pas_NO3/frac_pas' "${TMP}/tmp4.nc" "${OUT_DIR}/fertl_pasture_no3_${year}.nc"

  ${cdo} -O merge "${TMP}/manure_crops_${year}.nc" "${FRAC_CROP}" "${TMP}/tmp5.nc"
  ${cdo} -O expr,'manure_crops=Nmanure_Crop_App/frac_crop' "${TMP}/tmp5.nc" "${OUT_DIR}/manure_crops_${year}.nc"

  ${cdo} -O merge "${TMP}/manure_pasture_${year}.nc" "${FRAC_PAS}" "${TMP}/tmp6.nc"
  ${cdo} -O expr,'manure_pasture=Nmanure_Pas_App/frac_pas' "${TMP}/tmp6.nc" "${OUT_DIR}/manure_pasture_${year}.nc"

  ${cdo} -merge "${OUT_DIR}/fertl_crops_nh4_${year}.nc" \
  "${OUT_DIR}/fertl_crops_no3_${year}.nc" \
  "${OUT_DIR}/fertl_pasture_nh4_${year}.nc" \
  "${OUT_DIR}/fertl_pasture_no3_${year}.nc" \
  "${OUT_DIR}/manure_crops_${year}.nc" \
  "${OUT_DIR}/manure_pasture_${year}.nc" \
  "${OUT_DIR}/fert_man_icon_scaled_${year}.nc"

  ${cdo} -O expr,'fertl_crops_nh4=max(fertl_crops_nh4,0);fertl_crops_no3=max(fertl_crops_no3,0);fertl_pasture_nh4=max(fertl_pasture_nh4,0);fertl_pasture_no3=max(fertl_pasture_no3,0);manure_crops=max(manure_crops,0);manure_pasture=max(manure_pasture,0)' \
  "${OUT_DIR}/fert_man_icon_scaled_${year}.nc" \
  "${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"

done

for year in $(seq 1850 2024); do
  HARV_SRC="${LUH_LAND_USE_DIR}/bc_land_use_quincy_${year}.nc"
  HARV_TMP="${TMP}/harvest_fract_${year}.nc"
  FERT_ICON="${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"
  OUT_NMIP="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

  ${ncatted} -a _FillValue,fertl_crops_nh4,d,, -a missing_value,fertl_crops_nh4,d,, "${FERT_ICON}"
  ${ncatted} -a _FillValue,fertl_crops_no3,d,, -a missing_value,fertl_crops_no3,d,, "${FERT_ICON}"
  ${ncatted} -a _FillValue,fertl_pasture_nh4,d,, -a missing_value,fertl_pasture_nh4,d,, "${FERT_ICON}"
  ${ncatted} -a _FillValue,fertl_pasture_no3,d,, -a missing_value,fertl_pasture_no3,d,, "${FERT_ICON}"
  ${ncatted} -a _FillValue,manure_crops,d,, -a missing_value,manure_crops,d,, "${FERT_ICON}"
  ${ncatted} -a _FillValue,manure_pasture,d,, -a missing_value,manure_pasture,d,, "${FERT_ICON}"

  # harvest_fract z oryginalnego bc_land_use_quincy
  ${cdo} selvar,harvest_fract "${HARV_SRC}" "${HARV_TMP}"

  # merge: fertl_* + harvest_fract
  ${cdo} merge "${FERT_ICON}" "${HARV_TMP}" "${OUT_NMIP}"
done

for file in "${FINAL_DIR}"/bc_land_use_quincy_NMIP_????.nc; do

  [ -e "${file}" ] || continue

  echo "Reducing time dimension in ${file}"

  tmp_file="${file}.tmp.nc"

  # Dla jednego timestepu ncwa usunie wymiar time
  ${cdo} --reduce_dim copy "${file}" "${tmp_file}"

  # Podmiana pliku wynikowego
  mv "${tmp_file}" "${file}"
done

# create years 1700 till 1850 and 2025
# fertilization variables from 1850 year, harv from created earlier luh files 

# Exctract fert variables from year 1850
# Create 2025 year by copying fert data from 2024

for year in 1850 2024; do
  ${cdo} selvar,fertl_crops_nh4 "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/fertl_crops_nh4_${year}.nc
  ${cdo} selvar,fertl_crops_no3 "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/fertl_crops_no3_${year}.nc
  ${cdo} selvar,fertl_pasture_nh4 "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/fertl_pasture_nh4_${year}.nc
  ${cdo} selvar,fertl_pasture_no3 "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/fertl_pasture_no3_${year}.nc
  ${cdo} selvar,manure_crops    "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/manure_crops_${year}.nc
  ${cdo} selvar,manure_pasture  "${FINAL_DIR}"/bc_land_use_quincy_NMIP_${year}.nc "${TMP}"/manure_pasture_${year}.nc
done
# -----------------------------------------------------------
# Extract only harvest_fract for years 1700–1850
# ------------------------------------------------------------------
for year in $(seq 1700 1849); do
  HARV_SRC="${LUH_LAND_USE_DIR}/bc_land_use_quincy_${year}.nc"
  TMP_FILE="${TMP}/tmp_bc_land_use_quincy_NMIP_${year}.nc"
  OUT_NMIP="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

  # Select only harvest_fract variable
  ${cdo} selvar,harvest_fract "${HARV_SRC}" "${TMP_FILE}"
  ${cdo} merge "${TMP}"/fertl_crops_nh4_1850.nc "${TMP}"/fertl_crops_no3_1850.nc "${TMP}"/fertl_pasture_nh4_1850.nc \
  "${TMP}"/manure_crops_1850.nc \
  "${TMP}"/manure_pasture_1850.nc \
  "${TMP}"/fertl_pasture_no3_1850.nc \
  "${TMP_FILE}" "${OUT_NMIP}"
done

${cdo} selvar,harvest_fract "${LUH_LAND_USE_DIR}/bc_land_use_quincy_2025.nc" "${TMP}/harv_fract_2025.nc"
${cdo} merge "${TMP}/fertl_crops_nh4_2024.nc" "${TMP}/fertl_crops_no3_2024.nc" "${TMP}/fertl_pasture_nh4_2024.nc" "${TMP}/fertl_pasture_no3_2024.nc" "${TMP}/harv_fract_2025.nc" "${TMP}/manure_crops_2024.nc" "${TMP}/manure_pasture_2024.nc" "${FINAL_DIR}/bc_land_use_quincy_NMIP_2025.nc"

# An option to remove temporary directory
if [[ "$REMOVE_TMP" == "yes" ]]; then
  rm -f "${TMP}"
else
  echo "The "${TMP}" directory was not removed"
fi
