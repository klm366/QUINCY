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
#SBATCH --output=nfertil.o%j
#SBATCH --error=nfertil.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive
#SBATCH --ntasks=1
set -euo pipefail

module load cdo   # lub odpowiednik na Twoim klastrze
module load nco

FERT_DIR="/work/mj0143//icon_land_data_pool/trendy/2026/raw-data-20260702/nfert"
FRAC_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"
OUT_DIR="./FERT_final_3.8"
states="/work/mj0143/icon_land_data_pool/trendy/2026/raw-data-20260702/luh2/states_LUH2_GCB2026.nc"
GRID_ICON="/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc"
TMP="${OUT_DIR}"/tmp
#FINAL_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/luh/LUH_land_use_nmip_nfert_NEW"
LUH_LAND_USE_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/luh/LUH_land_use"
FINAL_DIR="${OUT_DIR}"/final
history="Calculate fertilisation application from NMIP"
mkdir -p "${OUT_DIR}"
mkdir -p "${TMP}"
mkdir -p "${FINAL_DIR}"

# Preparation of LUH file
# creating frac_crop and frac_pasture file !!! czeck with Chengs file later!!!!!!
cdo -O -selyear,1850/2024 ${states} ${TMP}/states_1850_2024.nc

# calculating crop and pasture from states
cdo -O -expr,'crop=c3ann+c4ann+c3per+c4per+c3nfx' \
  "${TMP}/states_1850_2024.nc" "${TMP}/crop.nc"
cdo -O -selname,pastr "${TMP}/states_1850_2024.nc" "${TMP}/pasture.nc"
#cdo -O -selname,range "${TMP}/states_1850_2025.nc" "${TMP}/range.nc"

cdo -O setmisstoc,0 "${TMP}/crop.nc" "${TMP}/crop0.nc"
cdo -O setmisstoc,0 "${TMP}/pasture.nc" "${TMP}/pasture0.nc"
#cdo -O setmisstoc,0 "${TMP}/range.nc" "${TMP}/range0.nc"

cdo -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/crop0.nc" "${TMP}/frac_crop.nc"
cdo -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/pasture0.nc" "${TMP}/frac_pas.nc"
#cdo -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/range0.nc" "${TMP}/frac_range.nc"

ncrename -O -h -v crop,frac_crop "${TMP}/frac_crop.nc"
ncrename -O -h -v pastr,frac_pas "${TMP}/frac_pas.nc"
#ncrename -O -h -v range,frac_range "${TMP}/frac_range.nc"

ncatted -O -h -a units,frac_crop,o,c,"1" "${TMP}/frac_crop.nc"
ncatted -O -h -a units,frac_pas,o,c,"1" "${TMP}/frac_pas.nc"
#ncatted -O -h -a units,frac_range,o,c,"1" "${TMP}/frac_range.nc"

#cdo -O merge "${TMP}/frac_crop.nc" "${TMP}/frac_pas.nc" "${TMP}/frac_range.nc" "${OUT_DIR}"/LUH2.nc 

###########################################################################


# Getting fertilization on the crop and pasture
cdo -mul "${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NH4.nc"
cdo -mul "${FERT_DIR}/Nfer_Crop_NO3.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NO3.nc"
cdo -mul "${FERT_DIR}/Nfer_Pas_NH4.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NH4.nc"
cdo -mul "${FERT_DIR}/Nfer_Pas_NO3.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NO3.nc"

# conversion from g N m-2 → kg N ha-1: *1e-3*1e4 = *10
# 
cdo -O -setmissval,nan "${TMP}/Nfer_Crop_NH4.nc" "${TMP}/Nfer_Crop_NH4_mvn.nc"
cdo -O -setmisstoc,0 "${TMP}/Nfer_Crop_NH4_mvn.nc" "${TMP}/Nfer_Crop_NH4_clean.nc"
cdo -O -mulc,10 "${TMP}/Nfer_Crop_NH4_clean.nc" "${TMP}/Nfer_Crop_NH4_kgNha.nc"

cdo -O -setmissval,nan "${TMP}/Nfer_Crop_NO3.nc" "${TMP}/Nfer_Crop_NO3_mvn.nc"
cdo -O -setmisstoc,0 "${TMP}/Nfer_Crop_NO3_mvn.nc" "${TMP}/Nfer_Crop_NO3_clean.nc"
cdo -O -mulc,10 "${TMP}/Nfer_Crop_NO3_clean.nc" "${TMP}/Nfer_Crop_NO3_kgNha.nc"

cdo -O -setmissval,nan "${TMP}/Nfer_Pas_NH4.nc" "${TMP}/Nfer_Pas_NH4_mvn.nc"
cdo -O -setmisstoc,0 "${TMP}/Nfer_Pas_NH4_mvn.nc" "${TMP}/Nfer_Pas_NH4_clean.nc"
cdo -O -mulc,10 "${TMP}/Nfer_Pas_NH4_clean.nc" "${TMP}/Nfer_Pas_NH4_kgNha.nc"

cdo -O -setmissval,nan "${TMP}/Nfer_Pas_NO3.nc" "${TMP}/Nfer_Pas_NO3_mvn.nc"
cdo -O -setmisstoc,0 "${TMP}/Nfer_Pas_NO3_mvn.nc" "${TMP}/Nfer_Pas_NO3_clean.nc"
cdo -O -mulc,10 "${TMP}/Nfer_Pas_NO3_clean.nc" "${TMP}/Nfer_Pas_NO3_kgNha.nc"

rm -f "${TMP}/Nfer_Crop_NH4_mvn.nc" "${TMP}/Nfer_Crop_NH4_clean.nc"
rm -f "${TMP}/Nfer_Crop_NO3_mvn.nc" "${TMP}/Nfer_Crop_NO3_clean.nc"
rm -f "${TMP}/Nfer_Pas_NH4_mvn.nc" "${TMP}/Nfer_Pas_NH4_clean.nc"
rm -f "${TMP}/Nfer_Pas_NO3_mvn.nc" "${TMP}/Nfer_Pas_NO3_clean.nc"

# Remappingto R2B5
cdo --no_history -remapcon,${GRID_ICON} "${TMP}/Nfer_Crop_NH4_kgNha.nc" "${TMP}/Nfer_Crop_NH4_kgNha_r2b5.nc"
cdo --no_history -remapcon,${GRID_ICON} "${TMP}/Nfer_Crop_NO3_kgNha.nc" "${TMP}/Nfer_Crop_NO3_kgNha_r2b5.nc"
cdo --no_history -remapcon,${GRID_ICON} "${TMP}/Nfer_Pas_NH4_kgNha.nc" "${TMP}/Nfer_Pas_NH4_kgNha_r2b5.nc"
cdo --no_history -remapcon,${GRID_ICON} "${TMP}/Nfer_Pas_NO3_kgNha.nc" "${TMP}/Nfer_Pas_NO3_kgNha_r2b5.nc"

# Splitting into years

cdo splityear "${TMP}/Nfer_Crop_NH4_kgNha_r2b5.nc" "${TMP}/Nfer_Crop_NH4_"
cdo splityear "${TMP}/Nfer_Crop_NO3_kgNha_r2b5.nc" "${TMP}/Nfer_Crop_NO3_"
cdo splityear "${TMP}/Nfer_Pas_NH4_kgNha_r2b5.nc" "${TMP}/Nfer_Pas_NH4_"
cdo splityear "${TMP}/Nfer_Pas_NO3_kgNha_r2b5.nc" "${TMP}/Nfer_Pas_NO3_"

# Going back to fert grid

for year in $(seq 1850 2024); do

  
  FRAC_FILE="${FRAC_DIR}/bc_land_frac_iq_13pfts_${year}.nc"
  FRAC_CROP="${TMP}/frac_crop_${year}.nc"
  FRAC_PAS="${TMP}/frac_pas_${year}.nc"

  # frac_crop = (pft11 + pft12) * fract_veg * land
  cdo expr,'frac_crop=(fract_pft11+fract_pft12)*fract_veg*fract_land*notsea' \
      "${FRAC_FILE}" "${FRAC_CROP}"

  # frac_pas = (pft09 + pft10) * fract_veg * land
  cdo expr,'frac_pas=(fract_pft09+fract_pft10)*fract_veg*fract_land*notsea' \
      "${FRAC_FILE}" "${FRAC_PAS}"

  # zabezpieczenia: frac==0 -> -1, a bardzo małe >0 -> 0.001
  cdo -O expr,'frac_crop=(frac_crop==0?-1:frac_crop)' \
    "${FRAC_CROP}" "${FRAC_CROP}.tmp"

  cdo -O expr,'frac_crop=(frac_crop>0 && frac_crop<=0.001?0.001:frac_crop)' \
    "${FRAC_CROP}.tmp" "${FRAC_CROP}"

  cdo -O expr,'frac_pas=(frac_pas==0?-1:frac_pas)' \
    "${FRAC_PAS}" "${FRAC_PAS}.tmp"

  cdo -O expr,'frac_pas=(frac_pas>0 && frac_pas<=0.001?0.001:frac_pas)' \
  "${FRAC_PAS}.tmp" "${FRAC_PAS}"

  cdo -O merge "${TMP}/Nfer_Crop_NH4_${year}.nc" "${FRAC_CROP}" "${TMP}/tmp1.nc"
  cdo -O expr,'fertl_crops_nh4=Nfer_Crop_NH4/frac_crop' "${TMP}/tmp1.nc" "${OUT_DIR}/fertl_crops_nh4_${year}.nc"

  cdo -O merge "${TMP}/Nfer_Crop_NO3_${year}.nc" "${FRAC_CROP}" "${TMP}/tmp2.nc"
  cdo -O expr,'fertl_crops_no3=Nfer_Crop_NO3/frac_crop' "${TMP}/tmp2.nc" "${OUT_DIR}/fertl_crops_no3_${year}.nc"

  cdo -O merge "${TMP}/Nfer_Pas_NH4_${year}.nc" "${FRAC_PAS}" "${TMP}/tmp3.nc"
  cdo -O expr,'fertl_pasture_nh4=Nfer_Pas_NH4/frac_pas' "${TMP}/tmp3.nc" "${OUT_DIR}/fertl_pasture_nh4_${year}.nc"

  cdo -O merge "${TMP}/Nfer_Pas_NO3_${year}.nc" "${FRAC_PAS}" "${TMP}/tmp4.nc"
  cdo -O expr,'fertl_pasture_no3=Nfer_Pas_NO3/frac_pas' "${TMP}/tmp4.nc" "${OUT_DIR}/fertl_pasture_no3_${year}.nc"

  cdo -merge "${OUT_DIR}/fertl_crops_nh4_${year}.nc" \
  "${OUT_DIR}/fertl_crops_no3_${year}.nc" \
  "${OUT_DIR}/fertl_pasture_nh4_${year}.nc" \
  "${OUT_DIR}/fertl_pasture_no3_${year}.nc" \
  "${OUT_DIR}/fert_man_icon_scaled_${year}.nc"

  cdo -O expr,'fertl_crops_nh4=max(fertl_crops_nh4,0);fertl_crops_no3=max(fertl_crops_no3,0);fertl_pasture_nh4=max(fertl_pasture_nh4,0);fertl_pasture_no3=max(fertl_pasture_no3,0)' \
  "${OUT_DIR}/fert_man_icon_scaled_${year}.nc" \
  "${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"

done

for year in $(seq 1850 2024); do
  HARV_SRC="${LUH_LAND_USE_DIR}/bc_land_use_quincy_${year}.nc"
  HARV_TMP="${TMP}/harvest_fract_${year}.nc"
  FERT_ICON="${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"
  OUT_NMIP="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

  ncatted -h -a _FillValue,fertl_crops_nh4,d,, -a missing_value,fertl_crops_nh4,d,, "${FERT_ICON}"
  ncatted -h -a _FillValue,fertl_crops_no3,d,, -a missing_value,fertl_crops_no3,d,, "${FERT_ICON}"
  ncatted -h -a _FillValue,fertl_pasture_nh4,d,, -a missing_value,fertl_pasture_nh4,d,, "${FERT_ICON}"
  ncatted -h -a _FillValue,fertl_pasture_no3,d,, -a missing_value,fertl_pasture_no3,d,, "${FERT_ICON}"

  # harvest_fract z oryginalnego bc_land_use_quincy
  cdo selvar,harvest_fract "${HARV_SRC}" "${HARV_TMP}"

  # merge: fertl_* + harvest_fract
  cdo merge "${FERT_ICON}" "${HARV_TMP}" "${OUT_NMIP}"
done

for file in "${FINAL_DIR}"/bc_land_use_quincy_NMIP_????.nc; do

  [ -e "${file}" ] || continue

  echo "Reducing time dimension in ${file}"

  tmp_file="${file}.tmp.nc"

  # Dla jednego timestepu ncwa usunie wymiar time
  cdo --reduce_dim copy "${file}" "${tmp_file}"

  # Podmiana pliku wynikowego
  mv "${tmp_file}" "${file}"
done
