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
cdo="cdo -s -O --no_history -b 64"

#----------------------------------------------------------------------

# path and name of the target output folder
FINAL_DIR=./test_nmips3
OUT_DIR="TEST_nmips3"
TMP="${OUT_DIR}"/tmp
REMOVE_TMP="yes"          # na czas testów

TEST_MODE="yes"          # yes = tylko TEST_YEARS, no = pełny przebieg
TEST_YEARS="1702 1850 2000 2025"

# lata z TEST_YEARS mieszczące się w zakresie [lo,hi] (albo pełny seq)
years_in_range() {
    local lo=$1 hi=$2 y
    if [[ "${TEST_MODE}" == "yes" ]]; then
        for y in ${TEST_YEARS}; do
            if (( y >= lo && y <= hi )); then echo "${y}"; fi
        done
    else
        seq "${lo}" "${hi}"
    fi
}

# rok nawozu: <1850 -> 1850, >2024 -> 2024
clamp() {
    local y=$1
    if (( y < 1850 )); then y=1850; fi
    if (( y > 2024 )); then y=2024; fi
    echo "${y}"
}

# lata, dla których trzeba policzyć fert_man_icon_scaled (1850 zawsze: BASE_1850)
if [[ "${TEST_MODE}" == "yes" ]]; then
    FERT_YEARS=$( { echo 1850; for y in ${TEST_YEARS}; do clamp "${y}"; done; } | sort -un )
else
    FERT_YEARS=$(seq 1850 2024)
fi
echo "FERT_YEARS: ${FERT_YEARS}"

SIMULATION_TYPE=NMIP
NMIP_RUN=NMIP_S2

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
${cdo} -O -selname,range "${TMP}/states_1850_2024.nc" "${TMP}/range.nc"

${cdo} -O setmisstoc,0 "${TMP}/crop.nc" "${TMP}/crop0.nc"
${cdo} -O setmisstoc,0 "${TMP}/pasture.nc" "${TMP}/pasture0.nc"
${cdo} -O setmisstoc,0 "${TMP}/range.nc" "${TMP}/range0.nc"

${cdo} -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/crop0.nc" "${TMP}/frac_crop.nc"
${cdo} -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/pasture0.nc" "${TMP}/frac_pas.nc"
${cdo} -O remapcon,"${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/range0.nc" "${TMP}/frac_range.nc"

ncrename -O -h -v crop,frac_crop "${TMP}/frac_crop.nc"
ncrename -O -h -v pastr,frac_pas "${TMP}/frac_pas.nc"
ncrename -O -h -v range,frac_range "${TMP}/frac_range.nc"

${ncatted} -O -a units,frac_crop,o,c,"1" "${TMP}/frac_crop.nc"
${ncatted} -O -a units,frac_pas,o,c,"1" "${TMP}/frac_pas.nc"
${ncatted} -O -h -a units,frac_range,o,c,"1" "${TMP}/frac_range.nc"
###########################################################################

# scale manure by crop and pasture fractions (frac_crop, frac_pas)
${cdo} -mul "${MANURE_DIR}/Nmanure_Crop_App.nc" "${TMP}/frac_crop.nc" "${TMP}/manure_crops_raw.nc"
${cdo} -mul "${MANURE_DIR}/Nmanure_Pas_App.nc" "${TMP}/frac_pas.nc" "${TMP}/manure_pasture_app_raw.nc"
${cdo} -mul "${MANURE_DIR}/Nmanure_Range_dep.nc" "${TMP}/frac_range.nc" "${TMP}/manure_range_dep_raw.nc"
${cdo} -mul "${MANURE_DIR}/Nmanure_Pas_dep.nc" "${TMP}/frac_pas.nc" "${TMP}/manure_pasture_dep_raw.nc"
# get fertilization on the crop and pasture
${cdo} -mul "${FERT_DIR}/Nfer_Crop_NH4.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NH4.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Crop_NO3.nc" "${TMP}/frac_crop.nc" "${TMP}/Nfer_Crop_NO3.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Pas_NH4.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NH4.nc"
${cdo} -mul "${FERT_DIR}/Nfer_Pas_NO3.nc" "${TMP}/frac_pas.nc" "${TMP}/Nfer_Pas_NO3.nc"

${cdo} -O -add "${TMP}/manure_pasture_dep_raw.nc" "${TMP}/manure_pasture_app_raw.nc" "${TMP}/manure_pasture_pas_raw.nc"
${cdo} -O -add "${TMP}/manure_pasture_pas_raw.nc" "${TMP}/manure_range_dep_raw.nc" "${TMP}/manure_pasture_raw.nc"
ncrename -O -h -v Nmanure_Pas_dep,manure_pasture_grid "${TMP}/manure_pasture_raw.nc"
# conversion from g N m-2 → kg N ha-1: *1e-3*1e4 = *10
# 
scale_to_kgNha() {
    local infile="$1"
    local outfile="$2"

    ${cdo} -O -setmissval,nan "${infile}" "${infile%.nc}_mvn.nc"
    ${cdo} -O -setmisstoc,0 "${infile%.nc}_mvn.nc" "${infile%.nc}_clean.nc"
    ${cdo} -O -mulc,10 "${infile%.nc}_clean.nc" "${outfile}"

    rm -f "${infile%.nc}_mvn.nc" "${infile%.nc}_clean.nc"
}

scale_to_kgNha "${TMP}/Nfer_Crop_NH4.nc" "${TMP}/Nfer_Crop_NH4_kgNha.nc"
scale_to_kgNha "${TMP}/Nfer_Crop_NO3.nc" "${TMP}/Nfer_Crop_NO3_kgNha.nc"
scale_to_kgNha "${TMP}/Nfer_Pas_NH4.nc"  "${TMP}/Nfer_Pas_NH4_kgNha.nc"
scale_to_kgNha "${TMP}/Nfer_Pas_NO3.nc"  "${TMP}/Nfer_Pas_NO3_kgNha.nc"
scale_to_kgNha "${TMP}/manure_crops_raw.nc"    "${TMP}/manure_crops_kgNha.nc"
scale_to_kgNha "${TMP}/manure_pasture_raw.nc"   "${TMP}/manure_pasture_kgNha.nc"

# remap to R2B5 ans plit years

remap_and_split() {
    local infile="$1"
    local prefix="$2"

    ${cdo} -remapcon,"${GRID_ICON}" "${infile}" "${infile%.nc}_r2b5.nc"
    ${cdo} splityear "${infile%.nc}_r2b5.nc" "${TMP}/${prefix}_"
}

remap_and_split "${TMP}/Nfer_Crop_NH4_kgNha.nc" "Nfer_Crop_NH4"
remap_and_split "${TMP}/Nfer_Crop_NO3_kgNha.nc" "Nfer_Crop_NO3"
remap_and_split "${TMP}/Nfer_Pas_NH4_kgNha.nc"  "Nfer_Pas_NH4"
remap_and_split "${TMP}/Nfer_Pas_NO3_kgNha.nc"  "Nfer_Pas_NO3"
remap_and_split "${TMP}/manure_crops_kgNha.nc"  "manure_crops"
remap_and_split "${TMP}/manure_pasture_kgNha.nc" "manure_pasture"



# going back to fert grid
for year in ${FERT_YEARS}; do
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
  ${cdo} -O expr,'manure_pasture=manure_pasture_grid/frac_pas' "${TMP}/tmp6.nc" "${OUT_DIR}/manure_pasture_${year}.nc"

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

FERT_VARS="fertl_crops_nh4,fertl_crops_no3,fertl_pasture_nh4,fertl_pasture_no3"
MANURE_VARS="manure_crops,manure_pasture"

BASE_1850="${OUT_DIR}/fert_man_icon_scaled_clean_1850.nc"

extract_harvest() {
    local year="$1"
    local outfile="$2"

    ${cdo} selvar,harvest_fract \
        "${LUH_LAND_USE_DIR}/bc_land_use_quincy_${year}.nc" \
        "${outfile}"
}

remove_fill_attributes() {
    local infile="$1"

    for var in fertl_crops_nh4 fertl_crops_no3 \
               fertl_pasture_nh4 fertl_pasture_no3 \
               manure_crops manure_pasture; do
        ${ncatted} -a _FillValue,"${var}",d,, \
                   -a missing_value,"${var}",d,, "${infile}"
    done
}
# frac crop/pasture for year L (same safety rules as in the main loop)
prep_frac() {
    local L="$1"
    local FRAC_FILE="${FRAC_DIR}/bc_land_frac_iq_13pfts_${L}.nc"

    ${cdo} expr,'frac_crop=(fract_pft11+fract_pft12)*fract_veg*fract_land*notsea' \
        "${FRAC_FILE}" "${TMP}/fc_${L}.nc"
    ${cdo} expr,'frac_pas=(fract_pft09+fract_pft10)*fract_veg*fract_land*notsea' \
        "${FRAC_FILE}" "${TMP}/fp_${L}.nc"

    # frac==0 -> -1 (later max(...,0) turns the result into 0); 0<frac<=0.001 -> 0.001
    ${cdo} expr,'frac_crop=(frac_crop==0?-1:((frac_crop>0 && frac_crop<=0.001)?0.001:frac_crop))' \
        "${TMP}/fc_${L}.nc" "${TMP}/frac_crop_L${L}.nc"
    ${cdo} expr,'frac_pas=(frac_pas==0?-1:((frac_pas>0 && frac_pas<=0.001)?0.001:frac_pas))' \
        "${TMP}/fp_${L}.nc" "${TMP}/frac_pas_L${L}.nc"
}

# data from year F divided by fractions from year L -> 6 variables in one file
scale_FL() {
    local F="$1" L="$2" out="$3"

    ${cdo} merge \
        "${TMP}/Nfer_Crop_NH4_${F}.nc" "${TMP}/Nfer_Crop_NO3_${F}.nc" \
        "${TMP}/Nfer_Pas_NH4_${F}.nc"  "${TMP}/Nfer_Pas_NO3_${F}.nc" \
        "${TMP}/manure_crops_${F}.nc"  "${TMP}/manure_pasture_${F}.nc" \
        "${TMP}/frac_crop_L${L}.nc"    "${TMP}/frac_pas_L${L}.nc" \
        "${TMP}/sfl.nc"

    ${cdo} expr,'fertl_crops_nh4=max(Nfer_Crop_NH4/frac_crop,0);fertl_crops_no3=max(Nfer_Crop_NO3/frac_crop,0);fertl_pasture_nh4=max(Nfer_Pas_NH4/frac_pas,0);fertl_pasture_no3=max(Nfer_Pas_NO3/frac_pas,0);manure_crops=max(Nmanure_Crop_App/frac_crop,0);manure_pasture=max(manure_pasture_grid/frac_pas,0)' \
        "${TMP}/sfl.nc" "${out}"
}
case "${NMIP_RUN}" in
    NMIP)
        # Fertilizer, manure and harvest all vary annually.
        for year in $(years_in_range 1850 2024); do
            BASE="${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"
            HARV_TMP="${TMP}/harvest_fract_${year}.nc"
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

            remove_fill_attributes "${BASE}"
            extract_harvest "${year}" "${HARV_TMP}"
            ${cdo} merge "${BASE}" "${HARV_TMP}" "${OUTFILE}"
        done
        ;;

    NMIP_S0)
        # Constant 1850 boundary condition; symlink it for 0850-2024.
        remove_fill_attributes "${BASE_1850}"
        extract_harvest 1850 "${TMP}/harvest_fract_1850.nc"
        ${cdo} merge "${BASE_1850}" "${TMP}/harvest_fract_1850.nc" \
            "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.nc"

        ${cdo} --reduce_dim copy \
            "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.nc" \
            "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.tmp.nc"
        mv "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.tmp.nc" \
           "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.nc"

        for year in $(years_in_range 850 2024); do
            if (( year == 1850 )); then continue; fi   # plik 1850 już istnieje
            printf -v year_padded '%04d' "${year}"
            ln -sf "${FINAL_DIR}/bc_land_use_quincy_NMIP_1850.nc" \
                "${FINAL_DIR}/bc_land_use_quincy_NMIP_${year_padded}.nc"
        done
        ;;

    NMIP_S2)
    # Fertilizer: year Y / frac Y.  Manure: year 1850 / frac Y (constant total per cell).
        for year in $(years_in_range 1850 2024); do
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

            prep_frac "${year}"
            scale_FL 1850 "${year}" "${TMP}/s2_man_full.nc"      # manure: F=1850, L=Y
            ${cdo} selvar,"${MANURE_VARS}" "${TMP}/s2_man_full.nc" "${TMP}/s2_manure.nc"

            BASE="${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"
            remove_fill_attributes "${BASE}"
            ${cdo} selvar,"${FERT_VARS}" "${BASE}" "${TMP}/s2_fert.nc"

            ${cdo} merge "${TMP}/s2_fert.nc" "${TMP}/s2_manure.nc" "${TMP}/s2_merged.nc"
            extract_harvest "${year}" "${TMP}/s2_harvest.nc"
            ${cdo} merge "${TMP}/s2_merged.nc" "${TMP}/s2_harvest.nc" "${OUTFILE}"
            remove_fill_attributes "${OUTFILE}"
    done
    ;;

    NMIP_S3)
        # Manure and harvest vary annually; fertilizer fixed at 1850.
        for year in $(years_in_range 1850 2024); do
            BASE="${OUT_DIR}/fert_man_icon_scaled_clean_${year}.nc"
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"
            remove_fill_attributes "${BASE}"
            remove_fill_attributes "${BASE_1850}"

            ${cdo} selvar,"${FERT_VARS}" "${BASE_1850}" "${TMP}/s3_fert.nc"
            ${cdo} selvar,"${MANURE_VARS}" "${BASE}" "${TMP}/s3_manure.nc"
            ${cdo} merge "${TMP}/s3_fert.nc" "${TMP}/s3_manure.nc" "${TMP}/s3_merged.nc"

            extract_harvest "${year}" "${TMP}/s3_harvest.nc"
            ${cdo} merge "${TMP}/s3_merged.nc" "${TMP}/s3_harvest.nc" "${OUTFILE}"
        done
        ;;

    NMIP_S6)
        # Fertilizer and manure: year Y / frac 1850 (land cover frozen at 1850).
        # Harvest fixed at 1850.
        extract_harvest 1850 "${TMP}/s6_harvest_1850.nc"
        prep_frac 1850

        for year in $(years_in_range 1850 2024); do
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

            scale_FL "${year}" 1850 "${TMP}/s6_fert_man.nc"   # F=Y, L=1850
            ${cdo} merge "${TMP}/s6_fert_man.nc" "${TMP}/s6_harvest_1850.nc" "${OUTFILE}"
            remove_fill_attributes "${OUTFILE}"
        done
        ;;
    NMIP_S11)
        # Fertilizer and manure fixed at 1850; harvest varies annually.
        remove_fill_attributes "${BASE_1850}"
        ${cdo} selvar,"${FERT_VARS},${MANURE_VARS}" "${BASE_1850}" "${TMP}/s11_base_1850.nc"

        for year in $(years_in_range 1850 2024); do
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

            extract_harvest "${year}" "${TMP}/s11_harvest.nc"
            ${cdo} merge "${TMP}/s11_base_1850.nc" "${TMP}/s11_harvest.nc" "${OUTFILE}"
        done
        ;;


    TRENDY)
        # fertilizer: clamp(year) do 1850-2024, harvest_fract: dany rok
        # (1700-1849 -> nawóz z 1850, 2025 -> nawóz z 2024)
        for fy in ${FERT_YEARS}; do
            remove_fill_attributes "${OUT_DIR}/fert_man_icon_scaled_clean_${fy}.nc"
        done

        for year in $(years_in_range 1700 2025); do
            fy=$(clamp "${year}")
            BASE="${OUT_DIR}/fert_man_icon_scaled_clean_${fy}.nc"
            OUTFILE="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"

            ${cdo} selvar,"${FERT_VARS}" "${BASE}" "${TMP}/trendy_fert_${year}.nc"
            extract_harvest "${year}" "${TMP}/trendy_harvest_${year}.nc"
            ${cdo} merge "${TMP}/trendy_fert_${year}.nc" \
                         "${TMP}/trendy_harvest_${year}.nc" "${OUTFILE}"
        done
        ;;
esac

#----------------------------------------------------------------------
# Remove the time dimension from regular annual output files
#----------------------------------------------------------------------
if [[ "${NMIP_RUN}" != "NMIP_S0" ]]; then
    for file in "${FINAL_DIR}"/bc_land_use_quincy_NMIP_????.nc; do
        [[ -e "${file}" ]] || continue

        echo "Reducing time dimension in ${file}"

        tmp_file="${file}.tmp.nc"
        ${cdo} --reduce_dim copy "${file}" "${tmp_file}"
        mv "${tmp_file}" "${file}"
    done
fi

#----------------------------------------------------------------------
# Optional cleanup
#----------------------------------------------------------------------
if [[ "${REMOVE_TMP}" == "yes" ]]; then
    rm -rf "${TMP}"
else
    echo "The ${TMP} directory was not removed"
fi