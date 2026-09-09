#!/bin/bash 


# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

grid="r2b5" # r2b5 or r2b4

# grid
if [[ "$grid" == "r2b5" ]]; then
  area="/work/mj0143/b301069/R2B5_grid_area.nc"
elif [[ "$grid" == "r2b4" ]]; then
  area="/work/mj0143/b301069/R2B4_grid_area.nc"
else
  echo "You have chosen grid does is not yet covered in the script"
fi

FINAL_DIR="." # directory where n depositions files are
outdir="./YOUROUTPUTFOLDER"
temp_path="${outdir}"/temp
FRAC_DIR="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"

mkdir -p "${outdir}"
mkdir -p "${temp_path}"
results_csv="${outdir}/fert_global_Tg.csv"
echo "year,var,value" > "${results_csv}"

for year in $(seq 1850 2024); do
	OUT_NMIP="${FINAL_DIR}/bc_land_use_quincy_NMIP_${year}.nc"
	FRAC_IQ="${FRAC_DIR}/bc_land_frac_iq_13pfts_${year}.nc"
	OUT_NMIP_kg="${temp_path}/bc_land_use_quincy_NMIP_${year}_kg.nc"
# 1) IQ frakcje crop/pasture
	${cdo} -O expr,'iq_fract_crop=(fract_pft11+fract_pft12)*fract_veg*land; iq_fract_past=(fract_pft09+fract_pft10)*fract_veg*land' \
  "${FRAC_IQ}" ${temp_path}/iq_frac.nc
	# cutting cell area
	${cdo} -O selvar,cell_area "${area}" ${temp_path}/cell_area.nc
	${cdo} -O merge "${OUT_NMIP}" ${temp_path}/iq_frac.nc \
	${temp_path}/cell_area.nc ${temp_path}/fert_new_merge.nc

	# kg N ha-1 → kg N m-2
  ${cdo} -O expr,'fert_crop_nh4_kg=(fertl_crops_nh4/10.0)*iq_fract_crop*cell_area;
	fert_crop_no3_kg=(fertl_crops_no3/10.0)*iq_fract_crop*cell_area;
	fert_past_nh4_kg=(fertl_pasture_nh4/10.0)*iq_fract_past*cell_area;
	fert_past_no3_kg=(fertl_pasture_no3/10.0)*iq_fract_past*cell_area' \
  ${temp_path}/fert_new_merge.nc "${OUT_NMIP_kg}"
	# Global sum for each var
	${cdo} -O fldsum -selname,fert_crop_nh4_kg,fert_crop_no3_kg,fert_past_nh4_kg,fert_past_no3_kg \
	"${OUT_NMIP_kg}" ${temp_path}/fertil_kg_global_${year}.nc
	# kg->Tg
	${cdo} -O mulc,1e-12 ${temp_path}/fertil_kg_global_${year}.nc ${temp_path}/fertil_Tg_global_${year}.nc
	echo "YEAR=${year}"
	tomp="${temp_path}/tmp_${year}.txt"
 ${cdo} outputtab,name,value ${temp_path}/fertil_Tg_global_${year}.nc > "${tomp}"
	while read var value; do
    echo "${year},${var},${value}" >> "${results_csv}"
  done < "${tomp}"
done
