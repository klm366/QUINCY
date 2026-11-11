#!/bin/bash 

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -O -s --no_history -b 64"

grid="r2b5" # r2b5 or r2b4
TEST_MODE="yes"      # yes = tylko TEST_YEARS, no = 1850-2024
TEST_YEARS="1702 1850 2000 2025"
NMIP_RUN_BUDGET="NMIP"   # TRENDY: bez manure; inne runy: z manure

if [[ "${TEST_MODE}" == "yes" ]]; then
    YEARS="${TEST_YEARS}"
else
    YEARS=$(seq 1850 2024)
fi
# grid
if [[ "$grid" == "r2b5" ]]; then
  area="/work/mj0143/b301069/R2B5_grid_area.nc"
elif [[ "$grid" == "r2b4" ]]; then
  area="/work/mj0143/b301069/R2B4_grid_area.nc"
else
  echo "You have chosen grid does is not yet covered in the script"
fi

#FINAL_DIR="/work/mj0143/b301069/data/NMIP3_input/IQ_grid/fert_man_new/S2" # directory where n depositions files are

FINAL_DIR="/work/mj0143/b301072//jsbach_repo/trendy/jsbach/scripts/preprocessing/fertilization//test_nmips1/" # directory where n depositions files are
outdir="./test-budget-nmips2"
temp_path="${outdir}"/temp
FRAC_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"
#FRAC_DIR="/work/mj0143/b301069/data/NMIP3_input/IQ_grid/LUH/output/"

mkdir -p "${outdir}"
mkdir -p "${temp_path}"
results_csv="${outdir}/fert_global_Tg.csv"
echo "year,var,value" > "${results_csv}"

for year in ${YEARS}; do
    OUT_NMIP="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"
    FRAC_IQ="${FRAC_DIR}/bc_land_frac_iq_13pfts_${year}.nc"
    OUT_NMIP_kg="${temp_path}/bc_land_use_quincy_NMIP_${year}_kg.nc"

    # skip years without input (e.g. no fraction file for 1702/2025)
    if [[ ! -e "${OUT_NMIP}" || ! -e "${FRAC_IQ}" ]]; then
        echo "SKIP ${year}: missing ${OUT_NMIP} or ${FRAC_IQ}"
        continue
    fi

    ${cdo} -O expr,'iq_fract_crop=(fract_pft11+fract_pft12)*fract_veg*land; iq_fract_past=(fract_pft09+fract_pft10)*fract_veg*land' \
        "${FRAC_IQ}" ${temp_path}/iq_frac.nc
    ${cdo} -O selvar,cell_area "${area}" ${temp_path}/cell_area.nc
    ${cdo} -O merge "${OUT_NMIP}" ${temp_path}/iq_frac.nc \
        ${temp_path}/cell_area.nc ${temp_path}/fert_new_merge.nc

    # fertilizer always; manure only if the variable is in the file
    EXPR='fert_crop_nh4_kg=(fertl_crops_nh4/10.0)*iq_fract_crop*cell_area;
          fert_crop_no3_kg=(fertl_crops_no3/10.0)*iq_fract_crop*cell_area;
          fert_past_nh4_kg=(fertl_pasture_nh4/10.0)*iq_fract_past*cell_area;
          fert_past_no3_kg=(fertl_pasture_no3/10.0)*iq_fract_past*cell_area'
    if ${cdo} showname "${OUT_NMIP}" | grep -qw manure_crops; then
        EXPR="${EXPR};
          manure_crops_kg=(manure_crops/10.0)*iq_fract_crop*cell_area;
          manure_past_kg=(manure_pasture/10.0)*iq_fract_past*cell_area"
    fi

    ${cdo} -O expr,"${EXPR}" ${temp_path}/fert_new_merge.nc "${OUT_NMIP_kg}"

    # global sum of all computed variables, kg -> Tg
    ${cdo} -O mulc,1e-12 -fldsum "${OUT_NMIP_kg}" ${temp_path}/fertil_Tg_global_${year}.nc

    echo "YEAR=${year}"
    tomp="${temp_path}/tmp_${year}.txt"
    ${cdo} outputtab,name,value ${temp_path}/fertil_Tg_global_${year}.nc > "${tomp}"
    # outputtab prints a header line starting with '#': skip it
    while read -r var value; do
        [[ "${var}" == \#* ]] && continue
        echo "${year},${var},${value}" >> "${results_csv}"
    done < "${tomp}"
done